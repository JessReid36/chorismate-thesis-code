#!/usr/bin/env python3
"""cs6_d48_prepare.py - DIAGNOSTIC 48: is the estimate S of DIAGNOSTIC 47 a rule? A test on more sites and on
arrangements, with basis functions at every positive charge.

Frame 42436, bare point charges, CAM-B3LYP-D3BJ/def2-TZVPD. Every run of parts A and B is a dianion that carries, at
each POSITIVE charge, the def2-TZVPD functions of a hydrogen atom and no atom (a ghost atom, as D47's level CAMTZG),
and is read against a charge-free control that carries the same ghost atoms.

  PART A  one positive charge at four outer sites: D47's ladder site, and the outer sites at which S of one charge
          turns negative at the smallest, the middle and the largest size. At each site, every size on a grid of
          STEP from the largest at which S says "does not take" at both states to the smallest at which S says
          "takes" at both. At R and TS.
  PART B  arrangements, at R and TS: the designs of D43 that D44 finds RELIABLE in def2-TZVPD; the design of D43
          whose positive charges S calls safe one by one and unsafe together; and, for each size of the ladder at
          which S of the single charge is below zero at both states, that charge with one equal negative charge
          at the outer site that raises S most.
  PART C  a calibration: one electron, the functions of the ghost atom and a bare charge of each size used, at the
          level of the runs without its dispersion term (CAM) and at Hartree-Fock (HF).
  PART V  the detachment energy of the pairs of part B, measured: the dianion and the monoanion with the pair's
          charges and no functions at them, from the committed inputs of D44, at R and TS.

S is formed here, before anything is run, from the committed maps of D47 by the formula of D47's pre-registration,
and is written into CONDITIONS.tsv with what it predicts.

select() and everything it uses need the standard library only, so the reader imports them and rebuilds every
table, input and point-charge file on its own.

Usage:  cs6_d48_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, re, shutil, sys

S6 = "phase2.2/cs_stage6"
GEOM = {"R": "d33a/sp_42436_R/substrate.xyz", "TS": "d33a/sp_42436_TS/substrate.xyz"}
SITES = "sites_f42436.tsv"
AMAT = "A_f42436_d33a.tsv"
FLOORS = "../cs_stage5/floors.tsv"
D33A_PTS = "d33a/points.xyz"
D47_PTS = "d47/points.xyz"
MAPS = {"di": "d47/T1_C_ctl_%s", "mo": "d47/T1_C_bare_%s", "di_sv": "d47/S1_C_ctl_%s", "mo_sv": "d47/S1_C_bare_%s"}
CTL, BARE = "d44/CAMTZ_C_ctl", "d44/CAMTZ_C_bare"
CTL_SV, BARE_SV = "d31/ctl", "d40/B_bare"
DI_INP = "d44/CAMTZ_B_n02_q0.50_di_R/job.inp"
MO_INP = "d44/CAMTZ_B_n02_q0.50_mo_R/job.inp"
D47_GCTL = "d47/CAMTZG_G_s1124_di_%s"
D47_P050 = "d47/CAMTZG_A_s1124_p050_di_%s"
NSITE, NPTS47 = 2626, 2629
STATES = ("R", "TS")
DI = "-2 1"
MO = "-1 2"
R_SAFE = 3.50
BOHR = 0.5291772083                        # A per bohr: the factor of cs5_amatrix.py and of D47
EV = 27.211386245988
RYD = EV / 2.0

LADDER_SITE = 1124
STEP = 5                                   # hundredths of e: the grid of sizes
S_BAND = 0.5                               # eV: S at or above +S_BAND predicts DID NOT TAKE, at or below -S_BAND TOOK ELECTRONS
R_MIN = 3.25                               # A: no two charges closer (G43.3)
VDE_MIN = 0.5                              # eV: the first-order detachment energy a pair must keep at both states
GHOST_OFF = 0.01                           # A: the ghost centre sits this far beyond the charge, on the line from the nearest atom
D44_DESIGNS = ["B_n02_q0.50", "B_n02_q0.75", "B_n02_q1.00", "B_n03_q0.75", "B_n05_q0.50", "B_n11_q0.25"]
RELIABLE44 = ["B_n02_q0.50", "B_n11_q0.25"]    # a check only: the rule below must give these
TOGETHER = ["B_n05_q0.50"]                     # a check only: the rule below must give this
SITES_A = [1124, 1166, 1025, 1407]             # a check only: the rules below must give these
D47_SIZES = ["0.75", "1.00"]                   # the sizes D47 ran with functions at the charge that D48 does not
ATOMS = ("C1", "H1", "H2", "C2", "C3", "O1", "O2", "O3", "C4", "H3", "C5", "H4", "C6", "C7", "H5", "C8", "H6", "C9", "H7",
         "O4", "H8", "C10", "O5", "O6")
NCOND = {"A": 17, "B": 5}                      # a check only

BODY = ["%maxcore 3000", "%pal nprocs 8 end", "%scf", "  MaxIter 300", "  STABPerform true", "end"]
PC_LINE = '%pointcharges "pc.pc"'
LINE = "! CAM-B3LYP D3BJ def2-TZVPD def2/J RIJCOSX TightSCF"
ONE = [("CAM", "! UKS CAM-B3LYP def2-TZVPD def2/J RIJCOSX TightSCF"), ("HF", "! UHF def2-TZVPD TightSCF")]
# Guess HCore: there is no atom for ORCA's default guess to take a density from; the core Hamiltonian holds the charge
ONE_BODY = ["%maxcore 3000", "%scf", "  MaxIter 300", "  Guess HCore", "end", PC_LINE, "* xyz -1 2", "H: %.6f 0.000000 0.000000" % GHOST_OFF, "*", ""]
ONE_TIMEOUT = 900                              # seconds: a calibration run is stopped after this
JOBS = [("A1", "cm_d48_a1"), ("A2", "cm_d48_a2"), ("A3", "cm_d48_a3"), ("A4", "cm_d48_a4"), ("B1", "cm_d48_b1"), ("V1", "cm_d48_v1"), ("E1", "cm_d48_e1")]
WITH_XYZ = ("ctl", "pdi", "pmo")                # the kinds of run whose input reads the geometry from substrate.xyz beside it
WITH_PC = ("cond", "one", "pdi", "pmo")         # the kinds of run that hold a point-charge file
WALL, MEM = "24:00:00", "48gb"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/d48"
OTHER = ["STAGE6_CRITERIA.txt", "STAGE6_AMENDMENT5_20261109.txt", "STAGE6_AMENDMENT6_20261117.txt", "STAGE6_AMENDMENT9_20261120.txt",
         "STAGE6_DIAGNOSTIC43_20261117.txt", "STAGE6_DIAGNOSTIC44_RESULT_20261120.txt",
         "STAGE6_DIAGNOSTIC47_20261121.txt", "STAGE6_DIAGNOSTIC47_RESULT_20261121.txt",
         "d43.tsv", "d44.tsv", "d47.tsv", "d43/CONDITIONS.tsv", "d47/CONDITIONS.tsv"]


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


def table(p):
    ls = [l for l in open(p).read().splitlines() if not l.startswith("#")]
    hdr = ls[0].split("\t")
    return [dict(zip(hdr, l.split("\t"))) for l in ls[1:] if l.strip()]


def read_geom(p):
    ls = open(p).read().split("\n")
    n = int(ls[0].split()[0])
    out = []
    for l in ls[2:2 + n]:
        f = l.split()
        out.append((f[0], float(f[1]), float(f[2]), float(f[3])))
    if len(out) != n:
        stop("%s does not hold the %d atoms its first line announces." % (p, n))
    return out


def read_col(p):
    return [float(l.split()[0]) for l in open(p).read().splitlines() if l.strip() and not l.startswith("#")]


def dist(p, q):
    return ((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2) ** 0.5


def pc_text(charges):
    """the committed format of d31, d33b and d43: a count, then ' %9.6f' and three 10-decimal coordinates per charge"""
    return "%d\n" % len(charges) + "".join(" %9.6f      %.10f      %.10f      %.10f\n" % (q, x[0], x[1], x[2]) for q, x in charges)


def pt_line(xyz):
    """one line of a points file for orca_vpot, in bohr, as d33a/points.xyz writes it"""
    return "%.10f %.10f %.10f" % (xyz[0] / BOHR, xyz[1] / BOHR, xyz[2] / BOHR)


def energy(p):
    """the last FINAL SINGLE POINT ENERGY of an output that ended normally"""
    t = open(p, errors="replace").read()
    e = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    if not e or "ORCA TERMINATED NORMALLY" not in t:
        stop("%s is not an output that ended normally with an energy." % p)
    return float(e[-1])


def vread(p, lines):
    """the potentials of one vpot.out, Eh per e, in point order; its points must be the lines given (bohr)"""
    L = [l.split() for l in open(p).read().splitlines() if l.strip()]
    if not L or len(L[0]) != 1 or int(L[0][0]) != len(lines) or len(L) != len(lines) + 1:
        stop("%s does not hold %d points." % (p, len(lines)))
    out = []
    for row, ref in zip(L[1:], lines):
        v = [float(t) for t in row]
        r = [float(t) for t in ref.split()]
        if len(v) != 4 or any(abs(v[j] - r[j]) > 1e-6 for j in range(3)):
            stop("the points of %s are not those of %s." % (p, D47_PTS))
        out.append(v[3])
    return out


def rname(cond, state):
    return "CAMTZG_%s_di_%s" % (cond, state)


def vname(cond, state, sp):
    """the run of part V: the pair's charges and no functions at them; sp is 'di' or 'mo'"""
    return "CAMTZ_%s_%s_%s" % (cond, sp, state)


def size_name(q):
    return "%03d" % round(100 * abs(q))


def s_terms(M, key, st, ch):
    """ch: [(q, xyz, row)]. Returns (first-order detachment energy eV, [S_i or None per charge]) by D47's formula."""
    mo, di = M[(key[1], st)], M[(key[0], st)]
    vde = M[("bare", key, st)] + sum(q * (mo[p] - di[p]) * EV for q, xyz, p in ch)
    out = []
    for i, (q, xyz, p) in enumerate(ch):
        if q <= 0:
            out.append(None)
            continue
        other = sum(qj * BOHR / dist(xyz, xj) for j, (qj, xj, pj) in enumerate(ch) if j != i)
        out.append(vde - EV * (mo[p] + other) - RYD * q * q)
    return vde, out


TZ, SV = ("di", "mo"), ("di_sv", "mo_sv")


def s_min(M, st, ch):
    vde, lst = s_terms(M, TZ, st, ch)
    k = min((s, i) for i, s in enumerate(lst) if s is not None)
    return k[0], k[1], vde


def predict(s):
    return "DID NOT TAKE" if s >= S_BAND else ("TOOK ELECTRONS" if s <= -S_BAND else "none")


def zero_size(M, st, p):
    """the size at which S of one positive charge alone at point p is zero"""
    v = (M[("mo", st)][p] - M[("di", st)][p]) * EV
    c = M[("bare", TZ, st)] - EV * M[("mo", st)][p]
    disc = v * v + 4 * RYD * c
    return None if disc < 0 else (v + disc ** 0.5) / (2 * RYD)


def select(s6):
    """Everything the diagnostic runs, from the committed files alone."""
    need = [GEOM["R"], GEOM["TS"], SITES, AMAT, FLOORS, D33A_PTS, D47_PTS, DI_INP, MO_INP] + OTHER
    for st in STATES:
        need += ["%s/vpot.out" % (MAPS[k] % st) for k in ("di", "mo", "di_sv", "mo_sv")]
        need += ["%s_%s/job.%s" % (c, st, e) for c in (CTL, BARE, CTL_SV, BARE_SV) for e in ("inp", "out")]
        need += ["%s/job.%s" % (D47_GCTL % st, e) for e in ("inp", "out")] + ["%s/job.out" % (D47_P050 % st)]
    for d in D44_DESIGNS:
        need += ["d43/%s_di_%s/pc.pc" % (d, st) for st in STATES]
    for f in need:
        if not os.path.isfile(os.path.join(s6, f)):
            stop("%s is missing." % f)
    # ---- the committed inputs are what this script expects, line for line
    want = dict(di=[LINE] + BODY + [PC_LINE, "* xyzfile %s substrate.xyz" % DI, ""], ctl=[LINE] + BODY + ["* xyzfile %s substrate.xyz" % DI, ""],
                mo=[LINE] + BODY + [PC_LINE, "* xyzfile %s substrate.xyz" % MO, ""])
    for k, f in (("di", DI_INP), ("ctl", "%s_R/job.inp" % CTL), ("mo", MO_INP)):
        if open(os.path.join(s6, f)).read().split("\n") != want[k]:
            stop("%s is not the committed input this script expects, line for line." % f)
    for st in STATES:
        if sha(os.path.join(s6, "%s_%s" % (CTL, st), "substrate.xyz")) != sha(os.path.join(s6, GEOM[st])) \
                or sha(os.path.join(s6, "%s_%s" % (CTL, st), "job.inp")) != sha(os.path.join(s6, "%s_R" % CTL, "job.inp")):
            stop("%s_%s is not a run of the committed %s geometry from the committed input." % (CTL, st, st))

    # ---- geometries, sites, coefficients (as D47 forms them)
    G = dict((st, read_geom(os.path.join(s6, GEOM[st]))) for st in STATES)
    for st in STATES:
        if [a[0] for a in G[st]] != [n.rstrip("0123456789") for n in ATOMS]:
            stop("the elements of %s are not those this script names, in order." % GEOM[st])
    fl = table(os.path.join(s6, FLOORS))
    if tuple(r.get("atom") for r in fl) != ATOMS:
        stop("floors.tsv does not name the 24 atoms as this script does.")
    rows = [l.split("\t") for l in open(os.path.join(s6, SITES)).read().splitlines() if not l.startswith("#")]
    head, rows = rows[0], rows[1:]
    c = dict((k, head.index(k)) for k in ("site", "layer", "x", "y", "z", "q_min_neg", "q_min_pos"))
    a = read_col(os.path.join(s6, AMAT))
    if len(rows) != NSITE or len(a) != NSITE:
        stop("the site table and the coefficient table do not both hold %d rows." % NSITE)
    plines = open(os.path.join(s6, D47_PTS)).read().split("\n")
    plines = [plines[0]] + [l for l in plines[1:] if l.strip()]
    if plines[0].strip() != str(NPTS47) or len(plines) != NPTS47 + 1:
        stop("%s does not hold the %d points of D47." % (D47_PTS, NPTS47))
    if plines[1:1 + NSITE] != [l for l in open(os.path.join(s6, D33A_PTS)).read().split("\n")[1:] if l.strip()]:
        stop("the first %d points of %s are not those of %s." % (NSITE, D47_PTS, D33A_PTS))

    def near(xyz):
        d = dict((st, min((dist(xyz, g[1:]), k) for k, g in enumerate(G[st]))) for st in STATES)
        m = min((d[st][0], d[st][1], st) for st in STATES)
        return dict(dR=d["R"][0], dTS=d["TS"][0], d=m[0], k=m[1], atom=ATOMS[m[1]], kst=m[2])

    sites = []
    for i, r in enumerate(rows):
        if int(r[c["site"]]) != i:
            stop("row %d of %s is not site %d." % (i, SITES, i))
        xyz = (float(r[c["x"]]), float(r[c["y"]]), float(r[c["z"]]))
        if pt_line(xyz) != plines[1 + i]:
            stop("site %d: its coordinates in bohr, formed here, are not line %d of %s." % (i, i + 2, D47_PTS))
        s = dict(i=i, xyz=xyz, qn=r[c["q_min_neg"]], qp=r[c["q_min_pos"]], a=a[i])
        s.update(near(xyz))
        sites.append(s)
    outer = [s for s in sites if s["dR"] >= R_SAFE and s["dTS"] >= R_SAFE and s["qn"] == "-1.00" and s["qp"] != "none"]
    lad = min(outer, key=lambda s: (s["a"], s["i"]))
    if lad["i"] != LADDER_SITE:
        stop("the outer site with the most negative coefficient is %d, not %d." % (lad["i"], LADDER_SITE))

    # ---- the committed maps of D47 and the committed energies under them
    M = {}
    for st in STATES:
        for k in ("di", "mo", "di_sv", "mo_sv"):
            M[(k, st)] = vread(os.path.join(s6, MAPS[k] % st, "vpot.out"), plines[1:])
        for key, cc, bb, run in ((TZ, CTL, BARE, ("di", "mo")), (SV, CTL_SV, BARE_SV, ("di_sv", "mo_sv"))):
            ec, eb = energy(os.path.join(s6, "%s_%s" % (cc, st), "job.out")), energy(os.path.join(s6, "%s_%s" % (bb, st), "job.out"))
            for kk, ref in zip(run, (ec, eb)):
                if abs(energy(os.path.join(s6, MAPS[kk] % st, "job.out")) - ref) > 1e-6:
                    stop("%s does not give back the committed energy of its control, so its map is not used." % (MAPS[kk] % st))
            M[("bare", key, st)] = (eb - ec) * EV

    def ch_of(pairs):
        return [(q, sites[i]["xyz"], i) for q, i in pairs]

    # ---- the estimate formed here is the one D47's reader committed
    n47 = 0
    for r in table(os.path.join(s6, "d47.tsv")):
        m = re.match(r"A_s(\d{4})_p(\d{3})$", r.get("condition", ""))
        if r.get("kind") == "cond" and m and r.get("S_CAMTZ_eV"):
            s = s_min(M, r["state"], ch_of([(int(m.group(2)) / 100.0, int(m.group(1)))]))[0]
            if abs(s - float(r["S_CAMTZ_eV"])) > 6e-5:
                stop("S of %s at %s formed here is %+.4f; d47.tsv holds %s." % (r["condition"], r["state"], s, r["S_CAMTZ_eV"]))
            n47 += 1
    if n47 < 20:
        stop("d47.tsv holds only %d single-charge rows with S to check this script's S against." % n47)

    conds = []

    def gpos(xyz):
        n = near(xyz)
        X = G[n["kst"]][n["k"]][1:]
        L = dist(X, xyz)
        return tuple(round(xyz[j] + GHOST_OFF * (xyz[j] - X[j]) / L, 6) for j in range(3))

    def add(name, part, pairs, where, gname, pc=None):
        if any(x["name"] == name for x in conds):
            stop("%s is selected twice." % name)
        ch = ch_of(pairs)
        charges = [(q, xyz) for q, xyz, p in ch]
        x = dict(name=name, part=part, pairs=list(pairs), charges=charges, rows=[p for q, xyz, p in ch], where=where, gctl=gname,
                 pc=pc if pc is not None else pc_text(charges), ghosts=[gpos(xyz) for q, xyz in charges if q > 0], S={}, Si={}, vde={}, names={}, pred={})
        for st in STATES:
            s, k, vde = s_min(M, st, ch)
            x["S"][st], x["names"][st], x["vde"][st], x["Si"][st] = s, pairs[k][1], vde, s_terms(M, TZ, st, ch)[1]
            x["pred"][st] = predict(s)
        for q, i in pairs:
            if sites[i] not in outer:
                stop("%s puts a charge at site %d, which is not an outer site." % (name, i))
        for u in range(len(charges)):
            for w in range(u + 1, len(charges)):
                if dist(charges[u][1], charges[w][1]) < R_MIN:
                    stop("%s puts two charges closer than %.2f A." % (name, R_MIN))
        conds.append(x)
        return x

    def both(pairs):
        return tuple(s_min(M, st, ch_of(pairs))[0] for st in STATES)

    # ---- PART A: the sites
    z = lambda s: min(zero_size(M, st, s["i"]) for st in STATES)
    order = sorted(outer, key=lambda s: (z(s), s["i"]))
    picks = [("the ladder site of D47", lad), ("smallest zero size", order[0]), ("middle zero size", order[(len(order) - 1) // 2]), ("largest zero size", order[-1])]
    if [s["i"] for w, s in picks] != SITES_A or len(set(SITES_A)) != 4:
        stop("the rules give the sites %s, not %s." % ([s["i"] for w, s in picks], SITES_A))
    asites = []
    for why, s in picks:
        grid = [k for k in range(STEP, 100 + STEP, STEP)]
        sv = dict((k, both([(k / 100.0, s["i"])])) for k in grid)
        lo = [k for k in grid if min(sv[k]) >= S_BAND]
        hi = [k for k in grid if max(sv[k]) <= -S_BAND]
        if not lo or not hi or max(lo) >= min(hi):
            stop("site %d has no size at which S says DID NOT TAKE below one at which it says TOOK ELECTRONS." % s["i"])
        sizes = [k for k in grid if max(lo) <= k <= min(hi)]
        if any(min(sv[k]) >= S_BAND for k in grid if k > max(lo)) or any(max(sv[k]) <= -S_BAND for k in grid if k < min(hi)):
            stop("S at site %d does not fall with the size as this script assumes." % s["i"])
        asites.append(dict(site=s["i"], why=why, sizes=sizes, zero=dict((st, zero_size(M, st, s["i"])) for st in STATES), atom=s["atom"], kst=s["kst"], d=s["d"], a=s["a"]))
        for k in sizes:
            add("A_s%04d_p%03d" % (s["i"], k), "A", [(k / 100.0, s["i"])], "site %d (%s)" % (s["i"], why), "G_s%04d" % s["i"])

    # ---- PART B: the committed designs
    d43c = dict((r["condition"], r) for r in table(os.path.join(s6, "d43/CONDITIONS.tsv")))
    d44t = table(os.path.join(s6, "d44.tsv"))
    des = {}
    for d in D44_DESIGNS:
        if d not in d43c or d43c[d].get("part") != "B":
            stop("%s is not a design of D43." % d)
        txt = open(os.path.join(s6, "d43/%s_di_R/pc.pc" % d)).read()
        if open(os.path.join(s6, "d43/%s_di_TS/pc.pc" % d)).read() != txt:
            stop("D43's two point-charge files of %s differ." % d)
        ls = [l.split() for l in txt.split("\n")[1:] if l.strip()]
        ids = [int(v) for v in d43c[d]["sites"].split(",")]
        if len(ids) != len(ls) or any(dist(sites[i]["xyz"], (float(f[1]), float(f[2]), float(f[3]))) > 1e-6 for i, f in zip(ids, ls)):
            stop("the charges of %s are not at the sites D43's table lists." % d)
        pairs = [(float(f[0]), i) for f, i in zip(ls, ids)]
        if pc_text([(q, sites[i]["xyz"]) for q, i in pairs]) != txt:
            stop("the point-charge file of %s is not what its sites and charges give." % d)
        w = [r.get("four_legs") for r in d44t if r.get("design") == d and r.get("level_key") == "CAMTZ"]
        if len(w) != 2:
            stop("d44.tsv does not hold the two rows of %s at CAMTZ." % d)
        des[d] = dict(pairs=pairs, txt=txt, reliable=(w == ["RELIABLE", "RELIABLE"]),
                      alone=all(min(both([(q, i)])) >= S_BAND for q, i in pairs if q > 0), together=max(both(pairs)) <= -S_BAND)
    rel = [d for d in D44_DESIGNS if des[d]["reliable"]]
    tog = [d for d in D44_DESIGNS if not des[d]["reliable"] and des[d]["alone"] and des[d]["together"]]
    if rel != RELIABLE44 or tog != TOGETHER:
        stop("the rules give the designs %s and %s, not %s and %s." % (rel, tog, RELIABLE44, TOGETHER))

    def gname_of(d, pairs):
        pos = [i for q, i in pairs if q > 0]
        return "G_s%04d" % pos[0] if len(pos) == 1 else "G_%s" % d[2:]

    for d, why in [(d, "D43 design, RELIABLE in D44 at def2-TZVPD") for d in rel] + [(d, "D43 design, each positive charge safe alone by S, unsafe together") for d in tog]:
        add(d, "B", des[d]["pairs"], "%s; sites %s" % (why, d43c[d]["sites"]), gname_of(d, des[d]["pairs"]), pc=des[d]["txt"])
    # ---- PART B: the ladder's charge with one equal negative charge
    pairsB = []
    for k in asites[0]["sizes"]:
        q = k / 100.0
        if max(both([(q, lad["i"])])) >= 0:
            continue
        best = None
        for s in outer:
            if s["i"] == lad["i"] or dist(s["xyz"], lad["xyz"]) < R_MIN:
                continue
            pr = [(q, lad["i"]), (-q, s["i"])]
            if min(s_min(M, st, ch_of(pr))[2] for st in STATES) < VDE_MIN:
                continue
            key = (-min(both(pr)), s["i"])
            if best is None or key < best[0]:
                best = (key, s, pr)
        if best is None:
            stop("no outer site takes the negative charge of the pair at %.2f e." % q)
        pairsB.append((k, best[1]["i"]))
        add("P_s%04d_p%03d_s%04d_m%03d" % (lad["i"], k, best[1]["i"], k), "B", best[2],
            "+%.2f e at site %d with -%.2f e at site %d, %.2f A away" % (q, lad["i"], q, best[1]["i"], dist(best[1]["xyz"], lad["xyz"])), "G_s%04d" % lad["i"])
    for part, n in NCOND.items():
        if sum(1 for x in conds if x["part"] == part) != n:
            stop("part %s has %d conditions, not the %d DIAGNOSTIC 48 was pre-registered with." % (part, sum(1 for x in conds if x["part"] == part), n))

    # ---- the ghost controls: one for each distinct set of ghost centres
    gctl = []
    for x in conds:
        old = [g for g in gctl if g["name"] == x["gctl"]]
        if old:
            if old[0]["ghosts"] != x["ghosts"]:
                stop("%s does not carry the ghost centres of its control %s." % (x["name"], x["gctl"]))
        else:
            gctl.append(dict(name=x["gctl"], ghosts=x["ghosts"], sites=[i for q, i in x["pairs"] if q > 0]))
    # ---- the sizes of the calibration: every size a positive charge has with functions at it, here or in D47
    sizes = sorted(set([round(100 * q) for x in conds for q, i in x["pairs"] if q > 0] + [round(100 * float(s)) for s in D47_SIZES]))
    sel = dict(conds=conds, gctl=gctl, G=G, sites=sites, outer=len(outer), asites=asites, pairsB=pairsB, sizes=sizes, M=M, lines=plines[1:],
               designs=dict((d, des[d]["pairs"]) for d in D44_DESIGNS))
    # ---- the ghost control of the ladder site is, line for line, the one D47 ran
    for st in STATES:
        if input_text(s6, sel, "gctl", None, st, "G_s%04d" % LADDER_SITE) != open(os.path.join(s6, D47_GCTL % st, "job.inp")).read():
            stop("the ghost control of site %d at %s is not the input D47 committed." % (LADDER_SITE, st))
    return sel


def input_text(s6, sel, kind, cond=None, state=None, gname=None, one=None):
    """the text of one job.inp. kind: 'ctl' (the plain committed control), 'gctl', 'cond', 'one', and of part V
    'pdi' and 'pmo' (the committed dianion and monoanion inputs of D44 with point charges, byte for byte)"""
    if kind == "one":
        return "\n".join([dict(ONE)[one]] + ONE_BODY)
    if kind == "ctl":
        return open(os.path.join(s6, "%s_R/job.inp" % CTL)).read()
    if kind in ("pdi", "pmo"):
        return open(os.path.join(s6, DI_INP if kind == "pdi" else MO_INP)).read()
    tz = open(os.path.join(s6, DI_INP if kind == "cond" else "%s_R/job.inp" % CTL)).read()
    ghosts = [g for g in sel["gctl"] if g["name"] == gname][0]["ghosts"] if kind == "gctl" else cond["ghosts"]
    ls = tz.split("\n")
    if ls[-2] != "* xyzfile %s substrate.xyz" % DI or ls[-1] != "":
        stop("the committed CAMTZ input does not end as this script expects.")
    geo = open(os.path.join(s6, GEOM[state])).read().split("\n")
    n = int(geo[0].split()[0])
    return "\n".join(ls[:-2] + ["* xyz %s" % DI] + geo[2:2 + n] + ["H: %.6f %.6f %.6f" % g for g in ghosts] + ["*", ""])


def one_pc(k):
    return pc_text([(k / 100.0, (0.0, 0.0, 0.0))])


def runs_of(sel):
    """every run in the order it is made within its job: dict(run, job, level, kind, cond, state, ref, gname, one, size)"""
    out = []

    def add(run, job, level, kind, cond=None, state="R", ref=None, gname=None, one=None, size=None):
        out.append(dict(run=run, job=job, level=level, kind=kind, cond=cond, state=state, ref=ref, gname=gname, one=one, size=size))

    cs = sel["conds"]
    cm = dict((x["name"], x) for x in cs)
    asites = sel["asites"]

    def block(job, gnames, names):
        add("%s_C_ctl_R" % job, job, "CAMTZ", "ctl")
        for g in gnames:
            for st in STATES:
                add("CAMTZG_%s_di_%s" % (g, st), job, "CAMTZG", "gctl", state=st, gname=g)
        for n in names:
            for st in STATES:
                add(rname(n, st), job, "CAMTZG", "cond", n, st, ref="CAMTZG_%s_di_%s" % (cm[n]["gctl"], st))

    single = lambda site: [x["name"] for x in cs if x["part"] == "A" and x["pairs"][0][1] == site]
    g0 = "G_s%04d" % asites[0]["site"]
    block("A1", [g0], single(asites[0]["site"]) + [x["name"] for x in cs if x["part"] == "B" and x["gctl"] == g0])
    for job, a in zip(("A2", "A3", "A4"), asites[1:]):
        block(job, ["G_s%04d" % a["site"]], single(a["site"]))
    rest = [x for x in cs if x["part"] == "B" and x["gctl"] != g0]
    gs = []
    for x in rest:
        if x["gctl"] not in gs:
            gs.append(x["gctl"])
    block("B1", gs, [x["name"] for x in sorted(rest, key=lambda x: (-min(x["S"].values()), x["name"]))])
    add("V1_C_ctl_R", "V1", "CAMTZ", "ctl")
    for k, j in sel["pairsB"]:
        n = "P_s%04d_p%03d_s%04d_m%03d" % (LADDER_SITE, k, j, k)
        if n not in cm:
            stop("the pair %s is not a condition." % n)
        for st in STATES:
            add(vname(n, st, "di"), "V1", "CAMTZ", "pdi", n, st)
            add(vname(n, st, "mo"), "V1", "CAMTZ", "pmo", n, st)
    for key, line in ONE:
        for k in sel["sizes"]:
            add("ONE_%s_p%03d" % (key, k), "E1", "ONE_%s" % key, "one", one=key, size=k)
    if len(set(r["run"] for r in out)) != len(out):
        stop("a run is listed twice.")
    if set(r["cond"] for r in out if r["kind"] == "cond") != set(cm) or set(r["gname"] for r in out if r["kind"] == "gctl") != set(g["name"] for g in sel["gctl"]):
        stop("the jobs do not hold every condition and every ghost control.")
    return out


def conditions_text(sel):
    L = ["condition\tpart\twhere\tsites\tcharges_e\tghost_control\tghost_centres\tS_R_eV\tS_TS_eV\tS_names_site_R\tS_names_site_TS"
         "\tfirst_order_VDE_R_eV\tfirst_order_VDE_TS_eV\tS_predicts_R\tS_predicts_TS\tpc_sha256"]
    for x in sel["conds"]:
        L.append("\t".join([x["name"], x["part"], x["where"], ",".join(str(i) for q, i in x["pairs"]), ",".join("%+.4f" % q for q, i in x["pairs"]),
                            x["gctl"], ";".join("%.6f %.6f %.6f" % g for g in x["ghosts"]),
                            "%+.4f" % x["S"]["R"], "%+.4f" % x["S"]["TS"], str(x["names"]["R"]), str(x["names"]["TS"]),
                            "%+.4f" % x["vde"]["R"], "%+.4f" % x["vde"]["TS"], x["pred"]["R"], x["pred"]["TS"],
                            hashlib.sha256(x["pc"].encode()).hexdigest()]))
    return "\n".join(L) + "\n"


def sites_text(sel):
    L = ["site\twhy\tnearest_atom\tnearest_state\tnearest_A\ta_kcal_per_e\tzero_size_R_e\tzero_size_TS_e\tsizes_e"]
    for a in sel["asites"]:
        L.append("\t".join([str(a["site"]), a["why"], a["atom"], a["kst"], "%.4f" % a["d"], "%+.6f" % a["a"], "%.4f" % a["zero"]["R"], "%.4f" % a["zero"]["TS"],
                            ",".join("%.2f" % (k / 100.0) for k in a["sizes"])]))
    return "\n".join(L) + "\n"


def runs_text(sel):
    L = ["run\tjob\tlevel\tkind\tcondition\tstate\treference_control\tsize_e"]
    for r in runs_of(sel):
        ref = r["ref"] or ("%s_R (committed)" % CTL if r["kind"] == "ctl" else ("%s_%s (committed)" % (CTL, r["state"]) if r["kind"] in ("gctl", "pdi") else ""))
        L.append("\t".join([r["run"], r["job"], r["level"], r["kind"], r["cond"] or (r["gname"] or ""), "" if r["kind"] == "one" else r["state"], ref,
                            "%.2f" % (r["size"] / 100.0) if r["size"] is not None else ""]))
    return "\n".join(L) + "\n"


def levels_text(sel):
    L = ["key\tlevel\tmethod_line\tinputs\tjobs"]
    L.append("\t".join(["CAMTZ", "CAM-B3LYP-D3BJ/def2-TZVPD", LINE,
                        "the committed control input, byte for byte; in V1 also the committed dianion and monoanion inputs of D44 with point charges, byte for byte",
                        "A1,A2,A3,A4,B1,V1 (the first run of each); V1 (the pairs without functions at the charges)"]))
    L.append("\t".join(["CAMTZG", "CAM-B3LYP-D3BJ/def2-TZVPD + a ghost H at every positive charge", LINE,
                        "the committed CAMTZ input with the committed geometry written out and one ghost atom H: for each positive charge after it", "A1,A2,A3,A4,B1"]))
    for key, line in ONE:
        L.append("\t".join(["ONE_%s" % key, "one electron, the def2-TZVPD functions of H as a ghost atom, a bare charge %.2f A from their centre" % GHOST_OFF, line,
                            "written by this script; no substrate, no atom, serial", "E1"]))
    return "\n".join(L) + "\n"


def sources(sel):
    out = [GEOM["R"], GEOM["TS"], SITES, AMAT, FLOORS, D33A_PTS, D47_PTS, DI_INP, MO_INP] + OTHER
    for st in STATES:
        for k in ("di", "mo", "di_sv", "mo_sv"):
            out += ["%s/vpot.out" % (MAPS[k] % st), "%s/job.out" % (MAPS[k] % st)]
        out += ["%s_%s/job.%s" % (c, st, e) for c in (CTL, BARE, CTL_SV, BARE_SV) for e in ("inp", "out")]
        out += ["%s/job.%s" % (D47_GCTL % st, e) for e in ("inp", "out")] + ["%s/job.out" % (D47_P050 % st)]
    for d in D44_DESIGNS:
        out += ["d43/%s_di_%s/pc.pc" % (d, st) for st in STATES]
    seen, uniq = set(), []
    for f in out:
        if f not in seen:
            seen.add(f)
            uniq.append(f)
    return uniq


def pbs_text(job, sel, runs, hpc):
    jk, pbs = job
    mine = [r for r in runs if r["job"] == jk]
    one = jk == "E1"
    L = ["#!/bin/bash", "#PBS -N %s" % pbs,
         "#PBS -l select=1:ncpus=8:mem=%s" % MEM, "#PBS -l walltime=%s" % WALL,
         "#PBS -m ae", "#PBS -M 18660916@sun.ac.za", "#PBS -j oe",
         "#PBS -o %s/run_%s.pbs.out" % (hpc, jk),
         "# Stage 6 DIAGNOSTIC 48: is the estimate S of DIAGNOSTIC 47 a rule?",
         "# This job is one of %d. It makes %d runs in the order below; a run that fails does not stop the ones after it." % (len(JOBS), len(mine)),
         "# No -q: the allocation is automatic and walltime %s routes to day. Restartable:" % WALL,
         "#   submitted again, it skips every run that ended normally.",
         'export PATH="/apps/openmpi/4.1.1/bin:$PATH"',
         'export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:'
         '/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"',
         "export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1",
         "ORCA=/home/apps2/ORCA/6.0.1",
         "W=%s" % hpc, 'cd "$W" || exit 1',
         'if ! mkdir .running_%s 2>/dev/null; then echo "FAIL: another job holds .running_%s in $W"; exit 1; fi' % (jk, jk),
         "trap 'rmdir \"$W/.running_%s\" 2>/dev/null' EXIT" % jk,
         'echo "host=$(hostname) start=$(date)"',
         'echo "free memory at start:"; free -g | head -2',
         "run() {",
         '  d="$W/$1"',
         '  if [ -s "$d/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$d/job.out"; then echo "$1 done before"',
         '  else ( cd "$d" && %s"$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); ' % ("timeout -k 10 %d " % ONE_TIMEOUT if one else "") +
         'echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi',
         "}"]
    for r in mine:
        L.append("run %s" % r["run"])
    L += ['echo "end=$(date)"',
          'n=0; for r in %s; do grep -q \'ORCA TERMINATED NORMALLY\' "$W/$r/job.out" 2>/dev/null && n=$((n+1)); done' % " ".join(r["run"] for r in mine),
          'echo "finished: $n of %d"' % len(mine)]
    return "\n".join(L) + "\n"


def nsum(sel, runs):
    n = 0
    for r in runs:
        n += 1                                                    # job.inp
        n += 1 if r["kind"] in WITH_XYZ else 0                    # substrate.xyz
        n += 1 if r["kind"] in WITH_PC else 0                     # pc.pc
    return n + 5 + len(JOBS)                                      # CONDITIONS, SITES_A, RUNS, LEVELS, SOURCES; the job scripts


def main():
    if len(sys.argv) != 4:
        stop("usage: cs6_d48_prepare.py <results-repo> <out-dir> <hpc-dir>")
    res, out, hpc = sys.argv[1:4]
    s6 = os.path.join(res, S6)
    if os.path.exists(out) and os.listdir(out):
        stop("%s exists and is not empty; this script writes into a new directory only." % out)
    sel = select(s6)
    runs = runs_of(sel)
    cmap = dict((x["name"], x) for x in sel["conds"])
    for f in sources(sel):
        if not os.path.isfile(os.path.join(s6, f)):
            stop("%s, which this diagnostic reads, is missing." % f)
    os.makedirs(out, exist_ok=True)
    sums = []
    for r in runs:
        d = os.path.join(out, r["run"])
        os.makedirs(d)
        x = cmap.get(r["cond"])
        with open(os.path.join(d, "job.inp"), "w") as f:
            f.write(input_text(s6, sel, r["kind"], x, r["state"], r["gname"], r["one"]))
        if r["kind"] in WITH_XYZ:
            shutil.copyfile(os.path.join(s6, GEOM[r["state"]]), os.path.join(d, "substrate.xyz"))
        if x is not None or r["kind"] == "one":
            with open(os.path.join(d, "pc.pc"), "w") as f:
                f.write(x["pc"] if x is not None else one_pc(r["size"]))
        for fn in sorted(os.listdir(d)):
            sums.append((sha(os.path.join(d, fn)), "%s/%s" % (r["run"], fn)))
    for fn, text in (("CONDITIONS.tsv", conditions_text(sel)), ("SITES_A.tsv", sites_text(sel)), ("RUNS.tsv", runs_text(sel)), ("LEVELS.tsv", levels_text(sel))):
        with open(os.path.join(out, fn), "w") as f:
            f.write(text)
    with open(os.path.join(out, "SOURCES.sha256"), "w") as f:
        for s in sources(sel):
            f.write("%s  %s\n" % (sha(os.path.join(s6, s)), s))
    for fn in ("CONDITIONS.tsv", "SITES_A.tsv", "RUNS.tsv", "LEVELS.tsv", "SOURCES.sha256"):
        sums.append((sha(os.path.join(out, fn)), fn))
    for job in JOBS:
        p = os.path.join(out, "run_%s.pbs" % job[0])
        with open(p, "w") as f:
            f.write(pbs_text(job, sel, runs, hpc))
        os.chmod(p, 0o755)
        sums.append((sha(p), "run_%s.pbs" % job[0]))
    if len(sums) != nsum(sel, runs):
        stop("%d files were written, not %d." % (len(sums), nsum(sel, runs)))
    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for h, n in sums:
            f.write("%s  %s\n" % (h, n))
    for x in sel["conds"]:
        print("   %-28s part %s  S R %+6.2f TS %+6.2f  %-14s %-14s  %s" % (x["name"], x["part"], x["S"]["R"], x["S"]["TS"], x["pred"]["R"], x["pred"]["TS"], x["where"]))
    for jk, pbs in JOBS:
        print("   job %s  %2d runs" % (jk, sum(1 for r in runs if r["job"] == jk)))
    print("   %d conditions, %d run directories, %d files checksummed, %d jobs" % (len(sel["conds"]), len(runs), len(sums), len(JOBS)))


if __name__ == "__main__":
    try:
        main()
    except Refused as e:
        sys.stderr.write("REFUSED  %s\n" % e)
        raise SystemExit(2)
