#!/usr/bin/env python3
"""cs6_d43_prepare.py - DIAGNOSTIC 43: bare charges at or beyond r_safe, on the committed grid.

74 single points on frame 42436 at CAM-B3LYP-D3BJ/def2-SVPD with the stability analysis.

  two controls   C_ctl_R, D31's committed charge-free dianion input at the reactant, byte for byte;
                 C_bare_R, D40's committed monoanion input for the bare substrate, byte for byte
  part B         six designs placed by the committed selection rule of D33.2 on the sites that lie
                 at or beyond 3.50 A from every substrate atom at both states, with one rule added:
                 at most half the charges of a design may be negative, so its charges sum to zero
                 or more
  part A         one bare charge of +1.00 e and one of -1.00 e at each of six such sites: the two
                 sites the designs lean on, and the four sites of D12's frozen sample that lie there

Every condition is run four times: the dianion and the monoanion, at the reactant and at the
transition state, with the same point charges. Every dianion input is the committed charged
input of D31, byte for byte; every monoanion input is that file with its charge-and-multiplicity
line changed from "-2 1" to "-1 2" and nothing else touched.

Nothing about a site or a design is typed in to be trusted. The distances are computed here from
the committed geometries. D12's sites are read from the table of its committed pre-registration.
The selection rule is first run on the whole grid without the added rule and must give back the
six committed designs of D33 stage B, site for site and charge for charge, with each score
within 1e-5 kcal/mol, before it is run on the outer sites. The lists typed in below are what
those rules give, and nothing is written unless they are exactly that.

select() and everything above it need the standard library only, so the reader imports them and
repeats the selection on its own. main() also loads numpy and the committed optimiser, and
refuses to go on unless the optimiser's own readers give the same sites, ranges, a and L as the
plain readers here, the same r_min, and the same designs when numpy orders the ranking.

With --plain that one cross-check is left out and nothing else changes: the same files are
written, byte for byte, on a machine that has no numpy or scipy. The committed inputs were
written without --plain.

Usage:  cs6_d43_prepare.py [--plain] <results-repo> <code-repo> <out-dir> <hpc-dir>
"""
import hashlib, math, os, re, shutil, sys

S6 = "phase2.2/cs_stage6"
GEOM = {"R": "d33a/sp_42436_R/substrate.xyz", "TS": "d33a/sp_42436_TS/substrate.xyz"}
TEMPLATE = "d31/B/n01_q0.25_R/job.inp"       # the committed charged dianion input; D31's gate 5 holds all 60 identical
CTL_DI = ("C_ctl_R", "d31/ctl_R")            # the charge-free dianion at the reactant
CTL_MO = ("C_bare_R", "d40/B_bare_R")        # the bare monoanion at the reactant
METHOD = "! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF"
DI, MONO = "-2 1", "-1 2"
STATES = ("R", "TS")
SPECIES = (("di", DI), ("mo", MONO))
R_SAFE = 3.50                                # A, D11
R_MIN = 3.25                                 # A, cs_opt_milp.py CONFIG["r_min"]; main() checks it against that file
CONV = 2.88                                  # e, the largest total D30 confirmed converged in def2-SVPD (D33.3)
NSITE = 2626
D12P = "STAGE6_DIAGNOSTIC12_20261027.txt"
D12_OUTER = [4, 396, 397, 1043]              # checked in select() against D12's table and the distances
ANCHOR_POS, ANCHOR_NEG = 1124, 557           # checked in select(): the outer both-sign sites of most negative and most positive a
LADDER = [(2, "0.50"), (2, "0.75"), (2, "1.00"), (3, "0.75"), (5, "0.50"), (11, "0.25")]
D33_DESIGNS = "designs_d33b.tsv"
SCORE_TOL = 1e-5                             # kcal/mol, on a committed score re-formed from the six-decimal matrix
SITES, AMAT, LMAT = "sites_f42436.tsv", "A_f42436_d33a.tsv", "L_f42436.tsv"
SOURCES = [SITES, AMAT, LMAT, D33_DESIGNS, D12P, GEOM["R"], GEOM["TS"], TEMPLATE,
           "d31/ctl_R/job.inp", "d31/ctl_R/job.out", "d31/ctl_TS/job.out",
           "d40/B_bare_R/job.inp", "d40/B_bare_R/job.out", "d40/B_bare_TS/job.out"]
NCOND, NRUN = 18, 74
WALL, MEM = "24:00:00", "48gb"


class Refused(Exception):
    pass


def stop(m):
    raise Refused(m)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def aname(site, q):
    """part A: one charge of +1.00 or -1.00 e at one site"""
    return "A_s%04d_%s100" % (site, "p" if q > 0 else "m")


def bname(n, q):
    """part B: a design of n charges of magnitude q; q is the text of the ladder"""
    return "B_n%02d_q%s" % (n, q)


def rname(cond, species, state):
    return "%s_%s_%s" % (cond, species, state)


def read_pc(p):
    """a committed-format point-charge file -> [(q, x, y, z)]"""
    ls = open(p).read().split("\n")
    n = int(ls[0].split()[0])
    out = [tuple(float(v) for v in l.split()) for l in ls[1:1 + n]]
    if len(out) != n or any(len(r) != 4 for r in out):
        stop("%s does not hold the %d charges its first line announces." % (p, n))
    return out


def read_geom(p):
    """an xyz file -> [(element, x, y, z)]"""
    ls = open(p).read().split("\n")
    n = int(ls[0].split()[0])
    out = []
    for l in ls[2:2 + n]:
        f = l.split()
        out.append((f[0], float(f[1]), float(f[2]), float(f[3])))
    if len(out) != n:
        stop("%s does not hold the %d atoms its first line announces." % (p, n))
    return out


def read_sites_plain(p):
    """the site table -> xyz, layer and admitted ranges, by the rule of cs_opt_milp.py read_sites"""
    rows = [l.rstrip("\n").split("\t") for l in open(p).read().splitlines() if not l.startswith("#")]
    head, rows = rows[0], rows[1:]
    c = dict((k, head.index(k)) for k in ("site", "layer", "x", "y", "z", "q_min_neg", "q_min_pos"))
    xyz, layer, ranges = [], [], []
    for k, r in enumerate(rows):
        if int(r[c["site"]]) != k:
            stop("row %d of %s is not site %d." % (k, p, k))
        xyz.append((float(r[c["x"]]), float(r[c["y"]]), float(r[c["z"]])))
        layer.append(r[c["layer"]])
        neg = r[c["q_min_neg"]] == "-1.00"
        pos = None if r[c["q_min_pos"]] == "none" else float(r[c["q_min_pos"]])
        if neg and pos is not None and pos == 0.0:
            rr = [(-1.0, 1.0)]
        else:
            rr = ([(-1.0, 0.0)] if neg else []) + ([(pos, 1.0)] if pos is not None else [])
        if not rr:
            stop("site %d admits no charge at all." % k)
        ranges.append(rr)
    return xyz, layer, ranges


def read_col(p):
    """a one-column matrix file -> its numbers, row i being site i"""
    return [float(l.split()[0]) for l in open(p).read().splitlines() if l.strip() and not l.startswith("#")]


def dist(p, q):
    return ((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2) ** 0.5


def nearest(pt, atoms):
    """(distance in A, element) of the atom nearest to a point"""
    best = (1e30, "?")
    for el, x, y, z in atoms:
        d = dist(pt, (x, y, z))
        if d < best[0]:
            best = (d, el)
    return best


def d12_sample(text):
    """the 36 sites of D12's frozen sample, from the table of its pre-registration"""
    out = []
    for l in text.splitlines():
        m = re.match(r"^\s+(\d+)\s+(cation|anion)\s+[A-Z][0-9A-Za-z]*\s+\d\.\d{3}\s+\d\.\d{3}\s+\d\.\d{3}\s+\d\.\d{3}\s+[A-Z][0-9A-Za-z]*\s+\d\.\d{3}", l)
        if m:
            out.append(int(m.group(1)))
    return out


def admits(ranges, i, v):
    return any(lo <= v <= hi for lo, hi in ranges[i])


def best_at(a, L, ranges, i, q):
    """D33.2: the better admissible sign of a charge of magnitude q at site i, by a q + L"""
    vals = dict((a[i] * v + L[i], v) for lo, hi in ranges[i] for v in (-q, q) if lo <= v <= hi)
    if not vals:
        return None, None
    m = min(vals)
    return m, vals[m]


def rule(xyz, a, L, ranges, pool, n, q, cap):
    """D33.2's greedy selection, as cs6_d33a.py holds it, on a pool of sites, passing over a
    negative charge once cap negative charges are taken. cap = n is D33.2 unchanged."""
    cand = []
    for i in pool:
        m, v = best_at(a, L, ranges, i, q)
        if m is not None:
            cand.append((m, i, v))
    cand.sort()
    sel, nneg = [], 0
    for m, i, v in cand:
        if len(sel) == n:
            break
        if v < 0 and nneg >= cap:
            continue
        if all(dist(xyz[i], xyz[j]) >= R_MIN for j, _, _ in sel):
            sel.append((i, v, m))
            nneg += 1 if v < 0 else 0
    if len(sel) != n:
        stop("the rule admits only %d sites for %d x %.2f." % (len(sel), n, q))
    return [s[0] for s in sel], [s[1] for s in sel], sum(s[2] for s in sel)


def select(s6):
    """Everything this diagnostic selects, from committed files alone. Standard library only."""
    for f in SOURCES:
        if not os.path.isfile(os.path.join(s6, f)):
            stop("%s is missing; the inputs are built from it." % f)
    xyz, layer, ranges = read_sites_plain(os.path.join(s6, SITES))
    a, L = read_col(os.path.join(s6, AMAT)), read_col(os.path.join(s6, LMAT))
    if len(xyz) != NSITE or len(a) != NSITE or len(L) != NSITE:
        stop("the grid is not %d sites with one a and one L each." % NSITE)
    head = " ".join(open(os.path.join(s6, AMAT)).read().split("\n")[:3])
    if sha(os.path.join(s6, SITES)) not in head:
        stop("%s was not built on the committed %s." % (AMAT, SITES))
    atoms = dict((st, read_geom(os.path.join(s6, GEOM[st]))) for st in STATES)
    if any(len(atoms[st]) != 24 for st in STATES):
        stop("a geometry does not hold the 24 atoms of the substrate.")
    near = dict((st, [nearest(xyz[i], atoms[st]) for i in range(NSITE)]) for st in STATES)
    dmin = [min(near["R"][i][0], near["TS"][i][0]) for i in range(NSITE)]
    outer = [i for i in range(NSITE) if dmin[i] >= R_SAFE]

    # the rule must give back D33 stage B's committed designs before it is used on anything else
    rows = [l.split("\t") for l in open(os.path.join(s6, D33_DESIGNS)).read().splitlines()
            if l.strip() and not l.startswith("#") and not l.startswith("name")]
    if len(rows) != 6:
        stop("%s does not list six designs." % D33_DESIGNS)
    for r in rows:
        sel, qs, sc = rule(xyz, a, L, ranges, range(NSITE), int(r[1]), float(r[2]), int(r[1]))
        # the committed scores were formed from the potentials before a was rounded to six decimals for its file
        if (",".join(str(s) for s in sel) != r[7] or ",".join("%+.4f" % c for c in qs) != r[8]
                or abs(sc - float(r[6])) > SCORE_TOL):
            stop("the selection rule written here does not give back the committed design %s." % r[0])

    # part A: the two sites the designs lean on, and D12's sites that lie at or beyond r_safe
    samp = d12_sample(open(os.path.join(s6, D12P)).read())
    if len(samp) != 36 or len(set(samp)) != 36 or 1690 not in samp:
        stop("the table of %s does not give D12's 36 sites." % D12P)
    d12o = sorted(i for i in samp if dmin[i] >= R_SAFE)
    if d12o != D12_OUTER:
        stop("D12's sites at or beyond %.2f A are %s, not %s." % (R_SAFE, d12o, D12_OUTER))
    both = [i for i in outer if admits(ranges, i, 1.0) and admits(ranges, i, -1.0)]
    apos = min(both, key=lambda i: (a[i], i))
    aneg = max(both, key=lambda i: (a[i], -i))
    if (apos, aneg) != (ANCHOR_POS, ANCHOR_NEG):
        stop("the outer both-sign sites of most negative and most positive a are %d and %d, not %d and %d."
             % (apos, aneg, ANCHOR_POS, ANCHOR_NEG))
    asites = [ANCHOR_POS, ANCHOR_NEG] + D12_OUTER
    if len(set(asites)) != 6 or any(i not in both for i in asites):
        stop("the six sites of part A are not six distinct outer sites that admit +1.00 e and -1.00 e.")

    cond = []          # (name, part, rule, sites, charges), in the order they are run
    for n, qt in LADDER:
        sel, qs, sc = rule(xyz, a, L, ranges, outer, n, float(qt), n // 2)
        cond.append((bname(n, qt), "B", "D33.2 on the outer sites, at most %d negative" % (n // 2), sel, qs))
    pairs = [c for c in cond if len(c[3]) == 2]
    if len(pairs) != 3 or any(c[3] != pairs[0][3] or [x > 0 for x in c[4]] != [x > 0 for x in pairs[0][4]] for c in pairs):
        stop("the three designs of two charges do not share their sites and signs.")
    if pairs[0][3] != [ANCHOR_NEG, ANCHOR_POS] or pairs[0][4][0] > 0 or pairs[0][4][1] < 0:
        stop("the design of two charges is not a negative charge at site %d and a positive one at site %d." % (ANCHOR_NEG, ANCHOR_POS))
    for s in asites:
        cond.append((aname(s, +1.0), "A", "D12" if s in D12_OUTER else "anchor", [s], [1.0]))
    for s in asites:
        cond.append((aname(s, -1.0), "A", "D12" if s in D12_OUTER else "anchor", [s], [-1.0]))
    if len(cond) != NCOND or len(set(c[0] for c in cond)) != NCOND:
        stop("%d conditions were listed, expected %d." % (len(cond), NCOND))
    for name, part, why, sel, qs in cond:
        if any(not admits(ranges, i, v) for i, v in zip(sel, qs)):
            stop("%s carries a charge its site does not admit." % name)
        if any(dmin[i] < R_SAFE for i in sel):
            stop("%s carries a charge inside %.2f A." % (name, R_SAFE))
        if part == "B":
            n = len(sel)
            if (sum(qs) < -5e-7 or sum(1 for v in qs if v < 0) > n // 2 or sum(abs(v) for v in qs) > CONV + 1e-9
                    or len(set(round(abs(v), 9) for v in qs)) != 1):
                stop("%s does not obey the rules of part B." % name)
            if any(dist(xyz[i], xyz[j]) < R_MIN for k, i in enumerate(sel) for j in sel[:k]):
                stop("%s holds two charges closer than r_min." % name)
    return dict(xyz=xyz, layer=layer, ranges=ranges, a=a, L=L, atoms=atoms, near=near, dmin=dmin, outer=outer,
                both=both, asites=asites, cond=cond)


def runs_of(cond):
    """the 74 runs in the order they are run: (run, condition or source, state, species, charge and multiplicity)"""
    runs = [(CTL_DI[0], CTL_DI[1], "R", None, None), (CTL_MO[0], CTL_MO[1], "R", None, None)]
    for c in cond:
        for sp_, cm in SPECIES:
            for st in STATES:
                runs.append((rname(c[0], sp_, st), c[0], st, sp_, cm))
    if len(runs) != NRUN or len(set(r[0] for r in runs)) != NRUN:
        stop("%d runs were listed, expected %d." % (len(runs), NRUN))
    return runs


def write_pc(path, pairs):
    """the committed format of d31 and d33b: a count, then ' %9.6f' and three 10-decimal coordinates"""
    L = ["%d" % len(pairs)]
    for q, (x, y, z) in pairs:
        L.append(" %9.6f      %.10f      %.10f      %.10f" % (q, x, y, z))
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")


def main():
    args = [x for x in sys.argv[1:] if x != "--plain"]
    plain = len(args) != len(sys.argv) - 1
    if len(args) != 4 or len(sys.argv) - 1 - len(args) > 1:
        stop("usage: cs6_d43_prepare.py [--plain] <results-repo> <code-repo> <out-dir> <hpc-dir>")
    res, code, out, hpc = args
    s6 = os.path.join(res, S6)
    if os.path.exists(out) and os.listdir(out):
        stop("%s exists and is not empty; this script writes into a new directory only." % out)
    sel = select(s6)
    xyz, layer, ranges, a, L = sel["xyz"], sel["layer"], sel["ranges"], sel["a"], sel["L"]
    atoms, near, dmin, cond = sel["atoms"], sel["near"], sel["dmin"], sel["cond"]

    # ---- the plain readers above against the committed optimiser's own, and r_min against its CONFIG
    if not plain:
        import importlib.util
        import numpy as np
        po = os.path.join(code, "phase2.2", "cs_opt_milp.py")
        if not os.path.isfile(po):
            stop("the committed optimiser is not at %s" % po)
        sp = importlib.util.spec_from_file_location("d43_opt", po)
        opt = importlib.util.module_from_spec(sp)
        sp.loader.exec_module(opt)
        oxyz, olayer, oranges = opt.read_sites(os.path.join(s6, SITES))
        oa = np.asarray(opt.read_matrix(os.path.join(s6, AMAT), NSITE, "A")).ravel()
        oL = np.asarray(opt.read_matrix(os.path.join(s6, LMAT), NSITE, "L")).ravel()
        if (not np.array_equal(np.asarray(oxyz), np.asarray(xyz)) or list(olayer) != layer
                or [[tuple(r) for r in rr] for rr in oranges] != ranges
                or not np.array_equal(oa, np.asarray(a)) or not np.array_equal(oL, np.asarray(L))):
            stop("the plain readers here do not give what the committed optimiser's readers give.")
        if abs(float(opt.CONFIG["r_min"]) - R_MIN) > 1e-12:
            stop("the committed r_min is %r, not %.2f." % (opt.CONFIG["r_min"], R_MIN))
        # numpy's ordering of the ranking, which is what cs6_d33a.py walked, must be the ordering used above
        for n, qt in LADDER:
            q = float(qt)
            vals = np.array([(best_at(a, L, ranges, i, q)[0] if (dmin[i] >= R_SAFE and best_at(a, L, ranges, i, q)[0] is not None) else 1e9)
                             for i in range(NSITE)])
            got, nneg = [], 0
            for i in np.argsort(vals):
                if vals[i] >= 1e9 or len(got) == n:
                    break
                v = best_at(a, L, ranges, int(i), q)[1]
                if v < 0 and nneg >= n // 2:
                    continue
                if all(np.linalg.norm(np.asarray(xyz[int(i)]) - np.asarray(xyz[j])) >= R_MIN for j in got):
                    got.append(int(i))
                    nneg += 1 if v < 0 else 0
            if got != [c for c in cond if c[0] == bname(n, qt)][0][3]:
                stop("numpy's ordering gives another design for %d x %s." % (n, qt))

    runs = runs_of(cond)
    for st in STATES:
        if sha(os.path.join(s6, "d31", "ctl_%s" % st, "substrate.xyz")) != sha(os.path.join(s6, GEOM[st])):
            stop("d31/ctl_%s/substrate.xyz is not frame 42436's committed %s geometry." % (st, st))
    tpl = open(os.path.join(s6, TEMPLATE)).read()
    tl = tpl.split("\n")
    if (tl[0] != METHOD or "STABPerform true" not in tpl or '%pointcharges "pc.pc"' not in tpl
            or sum(1 for l in tl if l.startswith("* xyzfile ")) != 1 or ("* xyzfile %s substrate.xyz" % DI) not in tl):
        stop("%s is not the charged CAM-B3LYP-D3BJ/def2-SVPD dianion input with the stability analysis." % TEMPLATE)
    for name, src in (CTL_DI, CTL_MO):
        t = open(os.path.join(s6, src, "job.inp")).read()
        want = DI if name == CTL_DI[0] else MONO
        if (t.split("\n")[0] != METHOD or "STABPerform true" not in t or "pointcharges" in t
                or ("* xyzfile %s substrate.xyz" % want) not in t.split("\n")):
            stop("%s/job.inp is not the charge-free input at charge and multiplicity %s." % (src, want))
        if sha(os.path.join(s6, src, "substrate.xyz")) != sha(os.path.join(s6, GEOM["R"])):
            stop("%s/substrate.xyz is not frame 42436's committed reactant geometry." % src)
        if "****ORCA TERMINATED NORMALLY****" not in open(os.path.join(s6, src, "job.out"), errors="replace").read():
            stop("%s/job.out did not terminate normally." % src)

    os.makedirs(out, exist_ok=True)
    sums = []
    cmap = dict((c[0], c) for c in cond)
    for run, src, st, sp_, cm in runs:
        d = os.path.join(out, run)
        os.makedirs(d)
        shutil.copyfile(os.path.join(s6, GEOM[st]), os.path.join(d, "substrate.xyz"))
        if sp_ is None:
            shutil.copyfile(os.path.join(s6, src, "job.inp"), os.path.join(d, "job.inp"))
        else:
            name, part, why, ss, qs = cmap[src]
            write_pc(os.path.join(d, "pc.pc"), [(v, xyz[i]) for i, v in zip(ss, qs)])
            new = [l.replace("* xyzfile %s " % DI, "* xyzfile %s " % cm) if l.startswith("* xyzfile ") else l for l in tl]
            with open(os.path.join(d, "job.inp"), "w") as f:
                f.write("\n".join(new))
        for fn in sorted(os.listdir(d)):
            sums.append((sha(os.path.join(d, fn)), "%s/%s" % (run, fn)))

    # ---- the checks, on what was actually written
    for run, src, st, sp_, cm in runs:
        d = os.path.join(out, run)
        if sha(os.path.join(d, "substrate.xyz")) != sha(os.path.join(s6, GEOM[st])):
            stop("%s/substrate.xyz is not frame 42436's committed %s geometry." % (run, st))
        if sp_ is None:
            if sorted(os.listdir(d)) != ["job.inp", "substrate.xyz"] or sha(os.path.join(d, "job.inp")) != sha(os.path.join(s6, src, "job.inp")):
                stop("%s is not the committed input of %s byte for byte." % (run, src))
            continue
        if sorted(os.listdir(d)) != ["job.inp", "pc.pc", "substrate.xyz"]:
            stop("%s holds files other than its input, its geometry and its point charges." % run)
        b = open(os.path.join(d, "job.inp")).read().split("\n")
        diff = [(x, y) for x, y in zip(tl, b) if x != y]
        if len(b) != len(tl):
            stop("%s: the input has a different number of lines from %s." % (run, TEMPLATE))
        if sp_ == "di":
            if diff or sha(os.path.join(d, "job.inp")) != sha(os.path.join(s6, TEMPLATE)):
                stop("%s is not %s byte for byte." % (run, TEMPLATE))
        elif len(diff) != 1 or diff[0] != ("* xyzfile %s substrate.xyz" % DI, "* xyzfile %s substrate.xyz" % MONO):
            stop("%s does not differ from %s in the charge-and-multiplicity line alone." % (run, TEMPLATE))
        name, part, why, ss, qs = cmap[src]
        pc = read_pc(os.path.join(d, "pc.pc"))
        if len(pc) != len(ss):
            stop("%s/pc.pc does not hold %d charges." % (run, len(ss)))
        for (q, x, y, z), i, v in zip(pc, ss, qs):
            if abs(q - v) > 1e-9 or max(abs(x - xyz[i][0]), abs(y - xyz[i][1]), abs(z - xyz[i][2])) > 1e-9:
                stop("%s/pc.pc does not hold site %d at %+.2f e." % (run, i, v))
            if min(nearest((x, y, z), atoms[s2])[0] for s2 in STATES) < R_SAFE:
                stop("%s/pc.pc holds a charge inside %.2f A." % (run, R_SAFE))
        if sha(os.path.join(d, "pc.pc")) != sha(os.path.join(out, rname(src, "di", "R"), "pc.pc")):
            stop("%s does not carry the same point-charge file as the other runs of %s." % (run, src))

    with open(os.path.join(out, "CONDITIONS.tsv"), "w") as f:
        f.write("condition\tpart\trule\tn\tq_e\tn_positive\tn_negative\tsum_q_e\tsum_abs_q_e\tclass\tsites\tcharges\t"
                "first_order_kcal\tL_kcal\tscore_kcal\tnearest_A_R\tnearest_A_TS\n")
        for name, part, why, ss, qs in cond:
            s = sum(qs)
            f.write("%s\t%s\t%s\t%d\t%.2f\t%d\t%d\t%+.4f\t%.4f\t%s\t%s\t%s\t%+.6f\t%+.6f\t%+.6f\t%.4f\t%.4f\n" % (
                name, part, why, len(ss), abs(qs[0]), sum(1 for v in qs if v > 0), sum(1 for v in qs if v < 0), s,
                sum(abs(v) for v in qs), "H" if s > -5e-7 else "N", ",".join(str(i) for i in ss),
                ",".join("%+.4f" % v for v in qs), math.fsum(a[i] * v for i, v in zip(ss, qs)), math.fsum(L[i] for i in ss),
                math.fsum([a[i] * v for i, v in zip(ss, qs)] + [L[i] for i in ss]),
                min(near["R"][i][0] for i in ss), min(near["TS"][i][0] for i in ss)))
    used = sorted(set(i for c in cond for i in c[3]))
    with open(os.path.join(out, "SITES.tsv"), "w") as f:
        f.write("site\tlayer\tx\ty\tz\ta_kcal_per_e\tL_kcal\tnearest_A_R\tnearest_atom_R\tnearest_A_TS\tnearest_atom_TS\tadmits\tused_in\n")
        for i in used:
            f.write("%d\t%s\t%.6f\t%.6f\t%.6f\t%+.6f\t%+.6f\t%.4f\t%s\t%.4f\t%s\t%s\t%s\n" % (
                i, layer[i], xyz[i][0], xyz[i][1], xyz[i][2], a[i], L[i], near["R"][i][0], near["R"][i][1],
                near["TS"][i][0], near["TS"][i][1], " ".join("[%.2f,%.2f]" % (lo, hi) for lo, hi in ranges[i]),
                ",".join(c[0] for c in cond if i in c[3])))
    with open(os.path.join(out, "SOURCES.sha256"), "w") as f:
        for s in SOURCES:
            f.write("%s  %s\n" % (sha(os.path.join(s6, s)), s))
    for fn in ("CONDITIONS.tsv", "SITES.tsv", "SOURCES.sha256"):
        sums.append((sha(os.path.join(out, fn)), fn))

    pbs = os.path.join(out, "run.pbs")
    with open(pbs, "w") as f:
        f.write("#!/bin/bash\n#PBS -N cm_d43\n")
        f.write("#PBS -l select=1:ncpus=8:mem=%s\n#PBS -l walltime=%s\n" % (MEM, WALL))
        f.write("#PBS -m ae\n#PBS -M 18660916@sun.ac.za\n#PBS -j oe\n")
        f.write("#PBS -o %s/run.pbs.out\n" % hpc)
        f.write("# Stage 6 DIAGNOSTIC 43: bare charges at or beyond r_safe, on the committed grid.\n")
        f.write("# C_ctl_R is D31's committed charge-free dianion input and C_bare_R is D40's committed\n")
        f.write("#   monoanion input for the bare substrate, each byte for byte: the two CONTROLS.\n")
        f.write("# Every other condition is run four times with one point-charge file: the dianion (di) and\n")
        f.write("#   the monoanion (mo), at the reactant and at the transition state.\n")
        f.write("# Order: the controls; the six designs of part B; the six charges of +1.00 e of part A;\n")
        f.write("#   then its six charges of -1.00 e. A run that fails does not stop the ones after it.\n")
        f.write("# No -q: the allocation is automatic and walltime %s routes to day. Restartable:\n" % WALL)
        f.write("#   submitted again, it skips every run that ended normally.\n")
        f.write('export PATH="/apps/openmpi/4.1.1/bin:$PATH"\n')
        f.write('export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:'
                '/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"\n')
        f.write("export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1\n")
        f.write("ORCA=/home/apps2/ORCA/6.0.1\n")
        f.write('W=%s\ncd "$W" || exit 1\n' % hpc)
        f.write('if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running in $W"; exit 1; fi\n')
        f.write("trap 'rmdir \"$W/.running\" 2>/dev/null' EXIT\n")
        f.write('echo "host=$(hostname) start=$(date)"\n')
        f.write('echo "free memory at start:"; free -g | head -2\n')
        f.write("run() {\n")
        f.write('  d="$W/$1"\n')
        f.write('  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"\n')
        f.write('  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); '
                'echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi\n')
        f.write("}\n")
        for r in runs:
            f.write("run %s\n" % r[0])
        f.write('echo "end=$(date)"\n')
        f.write('echo "finished: $(grep -l \'ORCA TERMINATED NORMALLY\' */job.out 2>/dev/null | wc -l) of %d"\n' % NRUN)
    os.chmod(pbs, 0o755)
    sums.append((sha(pbs), "run.pbs"))

    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for h, n in sums:
            f.write("%s  %s\n" % (h, n))
    print("   %d of %d sites lie at or beyond %.2f A from every atom at both states; the farthest is at %.3f A;"
          % (len(sel["outer"]), NSITE, R_SAFE, max(dmin)))
    print("   %d of them admit both +1.00 e and -1.00 e, and %d lie beyond 4.00 A"
          % (len(sel["both"]), sum(1 for i in sel["outer"] if dmin[i] > 4.0)))
    for name, part, why, ss, qs in cond:
        print("   %-14s %s  n %2d  sum %+.2f e  total %.2f e  first order %+8.4f  nearest %.3f A  sites %s" % (
            name, part, len(ss), sum(qs), sum(abs(v) for v in qs), math.fsum(a[i] * v for i, v in zip(ss, qs)),
            min(dmin[i] for i in ss), " ".join("%d%s" % (i, "+" if v > 0 else "-") for i, v in zip(ss, qs))))
    print("   %d run directories, %d files checksummed, %s, walltime %s" % (NRUN, len(sums), MEM, WALL))


if __name__ == "__main__":
    try:
        main()
    except Refused as e:
        sys.stderr.write("REFUSED  %s\n" % e)
        raise SystemExit(2)
