#!/usr/bin/env python3
"""cs6_d47_prepare.py - DIAGNOSTIC 47, a preliminary: when does a bare positive charge take electrons from the
substrate - its size, its place, its distance, and the reach of the basis set?

Frame 42436, bare point charges, dianion only. Four levels, the key being the name in the run directories:

  CAMSV    CAM-B3LYP-D3BJ/def2-SVPD     the working level
  CAMTZ    CAM-B3LYP-D3BJ/def2-TZVPD    the larger basis set of AMENDMENT 9
  CAMTZG   CAMTZ with the def2-TZVPD functions of a hydrogen atom, and no atom, placed at the charge (a ghost atom)
  CAMQZ    CAM-B3LYP-D3BJ/def2-QZVPPD   one basis step above CAMTZ

  PART A  the size ladder: one positive charge of each size of SIZES_A at site LADDER_SITE, at R and TS
  PART B  other surroundings: one positive charge of each size of SIZES_B at three sites picked by rule, at R
  PART C  the walk out: one charge of SIZE_C on the line from the nearest atom through LADDER_SITE, at the
          clearances of RUNGS, at R
  PART D1 reach, functions at the charge (CAMTZG): three of the conditions above again, with their own controls
  PART D2 reach, one basis step up (CAMQZ): the two designs of D43 that pass in D44's larger basis set, at R and TS
  PART E  the maps: orca_vpot at every committed site and every rung, on the charge-free dianion and monoanion
          controls of CAMSV and CAMTZ at R and TS

Every CAMSV and CAMTZ input is a committed input that has already run, byte for byte. A CAMQZ input is the
committed CAMTZ input with the basis set's name changed and nothing else. A CAMTZG input is the committed CAMTZ
input with its last line replaced by the committed geometry written out and one ghost-atom line after it.

select() and everything it uses need the standard library only, so the reader imports them and rebuilds every
table, input and point-charge file on its own.

Usage:  cs6_d47_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, os, shutil, sys

S6 = "phase2.2/cs_stage6"
GEOM = {"R": "d33a/sp_42436_R/substrate.xyz", "TS": "d33a/sp_42436_TS/substrate.xyz"}
SITES = "sites_f42436.tsv"
AMAT = "A_f42436_d33a.tsv"
FLOORS = "../cs_stage5/floors.tsv"
D33A_PTS = "d33a/points.xyz"
D33A_VPOT = {"R": "d33a/sp_42436_R/vpot.out", "TS": "d33a/sp_42436_TS/vpot.out"}
NSITE = 2626
STATES = ("R", "TS")
DI, MONO = "-2 1", "-1 2"
R_SAFE = 3.50
BOHR = 0.5291772083                        # A per bohr: the factor of cs5_amatrix.py (ORCA 6.0.1's own); it gives back every line of d33a/points.xyz

LADDER_SITE = 1124
SIZES_A = ["0.25", "0.50", "0.75", "1.00"]
SIZES_B = ["0.50", "1.00"]
CLASSES_B = [("carboxylate O", ("O1", "O2", "O5", "O6")), ("hydrogen", ("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8")),
             ("carbon", ("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"))]
SIZE_C = "1.00"
RUNGS = ["5.00", "7.00", "9.00"]           # A: clearance from every substrate atom at R and TS
GHOST_OFF = 0.01                           # A: the ghost centre sits this far beyond the charge, on the same line
GHOST = [("A", "0.50", ("R", "TS")), ("A", "0.75", ("R",)), ("C", "9.00", ("R",))]
DESIGNS_Q = ["B_n11_q0.25", "B_n02_q0.50"]
ATOMS = ("C1", "H1", "H2", "C2", "C3", "O1", "O2", "O3", "C4", "H3", "C5", "H4", "C6", "C7", "H5", "C8", "H6", "C9", "H7",
         "O4", "H8", "C10", "O5", "O6")
NCOND = {"A": 4, "B": 6, "C": 3, "Q": 2}   # a check only

BODY = ["%maxcore 3000", "%pal nprocs 8 end", "%scf", "  MaxIter 300", "  STABPerform true", "end"]
PC_LINE = '%pointcharges "pc.pc"'
LINE = {"CAMSV": "! CAM-B3LYP D3BJ def2-SVPD def2/J RIJCOSX TightSCF", "CAMTZ": "! CAM-B3LYP D3BJ def2-TZVPD def2/J RIJCOSX TightSCF",
        "CAMQZ": "! CAM-B3LYP D3BJ def2-QZVPPD def2/J RIJCOSX TightSCF"}
LINE["CAMTZG"] = LINE["CAMTZ"]
LEVELS = [
    dict(key="CAMSV", name="CAM-B3LYP-D3BJ/def2-SVPD", role="the working level",
         di="d31/B/n01_q0.25_R/job.inp", ctl="d31/ctl", bare="d40/B_bare"),
    dict(key="CAMTZ", name="CAM-B3LYP-D3BJ/def2-TZVPD", role="the larger basis set of AMENDMENT 9",
         di="d44/CAMTZ_B_n02_q0.50_di_R/job.inp", ctl="d44/CAMTZ_C_ctl", bare="d44/CAMTZ_C_bare"),
    dict(key="CAMTZG", name="CAM-B3LYP-D3BJ/def2-TZVPD + ghost H at the charge", role="reach: functions at the charge",
         di=None, ctl=None, bare=None),
    dict(key="CAMQZ", name="CAM-B3LYP-D3BJ/def2-QZVPPD", role="reach: one basis step up",
         di=None, ctl=None, bare=None),
]
JOBS = [("S1", "CAMSV", "cm_d47_s1"), ("T1", "CAMTZ", "cm_d47_t1"), ("T2", "CAMTZ", "cm_d47_t2"),
        ("G1", "CAMTZG", "cm_d47_g1"), ("Q1", "CAMQZ", "cm_d47_q1")]
WALL, MEM = "24:00:00", "48gb"
HPC = "/home/18660916/system_development/phase2.2/cs_stage6/d47"
OTHER = ["STAGE6_CRITERIA.txt", "STAGE6_AMENDMENT5_20261109.txt", "STAGE6_AMENDMENT9_20261120.txt",
         "STAGE6_DIAGNOSTIC43_20261117.txt", "STAGE6_DIAGNOSTIC43_RESULT_20261117.txt", "STAGE6_DIAGNOSTIC44_20261118.txt",
         "STAGE6_DIAGNOSTIC44_RESULT_20261120.txt", "d43.tsv", "d44.tsv", "d43/CONDITIONS.tsv", "d37/CAM_TZ_R/job.inp"]
D44_DESIGNS = ["B_n02_q0.50", "B_n02_q0.75", "B_n02_q1.00", "B_n03_q0.75", "B_n05_q0.50", "B_n11_q0.25"]


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


def rname(level, cond, state):
    return "%s_%s_di_%s" % (level, cond, state)


def select(s6):
    """Everything the diagnostic runs, from the committed files alone."""
    need = [GEOM["R"], GEOM["TS"], SITES, AMAT, FLOORS, D33A_PTS, D33A_VPOT["R"], D33A_VPOT["TS"]] + OTHER
    for lv in LEVELS[:2]:
        need += [lv["di"]] + ["%s_%s/job.%s" % (lv[c], st, e) for c in ("ctl", "bare") for st in STATES for e in ("inp", "out")]
    for d in DESIGNS_Q:
        need += ["d43/%s_di_%s/%s" % (d, st, f) for st in STATES for f in ("pc.pc", "job.out")]
    for d in D44_DESIGNS:
        need += ["d44/CAMTZ_%s_di_%s/%s" % (d, st, f) for st in STATES for f in ("pc.pc", "job.out")]
    for f in need:
        if not os.path.isfile(os.path.join(s6, f)):
            stop("%s is missing." % f)
    # ---- the committed inputs are what this script expects, line for line
    for lv in LEVELS[:2]:
        want = dict(di=[LINE[lv["key"]]] + BODY + [PC_LINE, "* xyzfile %s substrate.xyz" % DI, ""],
                    ctl=[LINE[lv["key"]]] + BODY + ["* xyzfile %s substrate.xyz" % DI, ""],
                    bare=[LINE[lv["key"]]] + BODY + ["* xyzfile %s substrate.xyz" % MONO, ""])
        got = dict(di=lv["di"], ctl="%s_R/job.inp" % lv["ctl"], bare="%s_R/job.inp" % lv["bare"])
        for k in want:
            if open(os.path.join(s6, got[k])).read().split("\n") != want[k]:
                stop("%s is not the committed input this script expects, line for line." % got[k])
        for c in ("ctl", "bare"):
            for st in STATES:
                if sha(os.path.join(s6, "%s_%s" % (lv[c], st), "substrate.xyz")) != sha(os.path.join(s6, GEOM[st])) \
                        or sha(os.path.join(s6, "%s_%s" % (lv[c], st), "job.inp")) != sha(os.path.join(s6, "%s_R" % lv[c], "job.inp")):
                    stop("%s_%s is not a run of the committed %s geometry from the committed input." % (lv[c], st, st))
    if open(os.path.join(s6, "d37/CAM_TZ_R/job.inp")).read().split("\n")[0] != LINE["CAMTZ"]:
        stop("line 1 of d37/CAM_TZ_R/job.inp is not the method line of CAMTZ.")
    if LINE["CAMQZ"] != LINE["CAMTZ"].replace("def2-TZVPD", "def2-QZVPPD") or LINE["CAMQZ"] == LINE["CAMTZ"]:
        stop("the method line of CAMQZ is not that of CAMTZ with the basis set's name changed.")

    # ---- geometries, sites, coefficients
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
    plines = open(os.path.join(s6, D33A_PTS)).read().split("\n")
    if plines[0].strip() != str(NSITE) or len([l for l in plines[1:] if l.strip()]) != NSITE:
        stop("%s does not hold the %d sites." % (D33A_PTS, NSITE))

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
            stop("site %d: its coordinates in bohr, formed here, are not line %d of %s." % (i, i + 2, D33A_PTS))
        s = dict(i=i, xyz=xyz, qn=r[c["q_min_neg"]], qp=r[c["q_min_pos"]], a=a[i])
        s.update(near(xyz))
        sites.append(s)
    outer = [s for s in sites if s["dR"] >= R_SAFE and s["dTS"] >= R_SAFE and s["qn"] == "-1.00" and s["qp"] != "none"]
    lad = min(outer, key=lambda s: (s["a"], s["i"]))
    if lad["i"] != LADDER_SITE:
        stop("the outer site with the most negative coefficient is %d, not %d." % (lad["i"], LADDER_SITE))

    conds = []

    def add(name, part, charges, where, info, levels):
        if any(x["name"] == name for x in conds):
            stop("%s is selected twice." % name)
        x = dict(name=name, part=part, charges=charges, pc=pc_text(charges), where=where, levels=levels, ghost=None)
        x.update(info)
        conds.append(x)

    def info_site(s, q):
        return dict(site=s["i"], q=q, xyz=s["xyz"], a=s["a"], first=s["a"] * q, atom=s["atom"], kst=s["kst"], d=s["d"], pt=s["i"])

    # ---- PART A
    for size in SIZES_A:
        q = float(size)
        add("A_s%04d_p%03d" % (lad["i"], round(100 * q)), "A", [(q, lad["xyz"])], "site %d" % lad["i"], info_site(lad, q),
            {"CAMSV": STATES, "CAMTZ": STATES})
    # ---- PART B
    bsites = []
    for cname, atoms in CLASSES_B:
        cand = [s for s in outer if s["atom"] in atoms and s["i"] != lad["i"]]
        if not cand:
            stop("no outer site has a %s as its nearest atom." % cname)
        b = min(cand, key=lambda s: (s["a"], s["i"]))
        bsites.append((cname, b, len(cand)))
        for size in SIZES_B:
            q = float(size)
            add("A_s%04d_p%03d" % (b["i"], round(100 * q)), "B", [(q, b["xyz"])], "site %d (%s)" % (b["i"], cname), info_site(b, q),
                {"CAMSV": ("R",), "CAMTZ": ("R",)})
    # ---- PART C: the line of D11, from the nearest atom through the site, at fixed clearances
    X = G[lad["kst"]][lad["k"]][1:]
    L = dist(X, lad["xyz"])
    u = tuple((lad["xyz"][j] - X[j]) / L for j in range(3))
    rungs = []
    for r in RUNGS:
        want, t = float(r), L
        while near(tuple(X[j] + t * u[j] for j in range(3)))["d"] < want:
            t += 0.0005
            if t > 40:
                stop("the line does not reach a clearance of %s A." % r)
        lo, hi = t - 0.0005, t
        for _ in range(60):
            m = 0.5 * (lo + hi)
            if near(tuple(X[j] + m * u[j] for j in range(3)))["d"] >= want:
                hi = m
            else:
                lo = m
        p = tuple(round(X[j] + hi * u[j], 6) for j in range(3))
        n = near(p)
        if abs(n["d"] - want) > 2e-4:
            stop("the rung at %s A has a clearance of %.4f A." % (r, n["d"]))
        q = float(SIZE_C)
        rungs.append((r, p))
        add("L_s%04d_r%03d_p%03d" % (lad["i"], round(100 * want), round(100 * q)), "C", [(q, p)], "line of site %d, clearance %s A" % (lad["i"], r),
            dict(site=None, q=q, xyz=p, a=None, first=None, atom=n["atom"], kst=n["kst"], d=n["d"], pt=NSITE + len(rungs) - 1),
            {"CAMSV": ("R",), "CAMTZ": ("R",)})
    # ---- PART D1: the ghost centre, 0.01 A beyond the charge on the line
    gctl = []
    for part, tag, states in GHOST:
        if part == "A":
            x = [y for y in conds if y["part"] == "A" and abs(y["q"] - float(tag)) < 1e-9][0]
            gk = "G_s%04d" % lad["i"]
        else:
            x = [y for y in conds if y["part"] == "C" and y["where"].endswith("clearance %s A" % tag)][0]
            gk = "G_r%03d" % round(100 * float(tag))
        gpos = tuple(round(x["xyz"][j] + GHOST_OFF * u[j], 6) for j in range(3))
        x["ghost"], x["gctl"] = gpos, gk
        x["levels"] = dict(x["levels"])
        x["levels"]["CAMTZG"] = states
        old = [g for g in gctl if g["name"] == gk]
        if old:
            old[0]["states"] = tuple(st for st in STATES if st in old[0]["states"] or st in states)
        else:
            gctl.append(dict(name=gk, ghost=gpos, states=tuple(states), pt=x["pt"]))
    # ---- PART D2: two committed designs
    d43c = dict((r["condition"], r) for r in table(os.path.join(s6, "d43/CONDITIONS.tsv")))
    for d in DESIGNS_Q:
        if d not in d43c or d43c[d].get("part") != "B":
            stop("%s is not a design of D43." % d)
        txt = open(os.path.join(s6, "d43/%s_di_R/pc.pc" % d)).read()
        if open(os.path.join(s6, "d43/%s_di_TS/pc.pc" % d)).read() != txt:
            stop("D43's two point-charge files of %s differ." % d)
        ls = [l.split() for l in txt.split("\n")[1:] if l.strip()]
        ch = [(float(f[0]), (float(f[1]), float(f[2]), float(f[3]))) for f in ls]
        ids = [int(v) for v in d43c[d]["sites"].split(",")]
        if len(ids) != len(ch) or any(dist(sites[i]["xyz"], xq[1]) > 1e-6 for i, xq in zip(ids, ch)):
            stop("the charges of %s are not at the sites D43's table lists." % d)
        add(d, "Q", ch, "D43 design, sites %s" % d43c[d]["sites"], dict(site=None, q=None, xyz=None, a=None, first=None, atom=None, kst=None, d=None, pt=None, ids=ids),
            {"CAMQZ": STATES})
        conds[-1]["pc"] = txt
    for part, n in NCOND.items():
        if sum(1 for x in conds if x["part"] == part) != n:
            stop("part %s does not have the %d conditions DIAGNOSTIC 47 was pre-registered with." % (part, n))
    points = [s["xyz"] for s in sites] + [p for r, p in rungs]
    return dict(conds=conds, gctl=gctl, G=G, sites=sites, outer=len(outer), bsites=[(cn, b["i"], n) for cn, b, n in bsites], rungs=rungs,
                points=points, line=dict(atom=lad["atom"], state=lad["kst"], u=u, L=L), levels=[dict(l) for l in LEVELS])


def input_text(s6, sel, level, kind, cond=None, state=None, gname=None):
    """the text of one job.inp. kind: 'cond', 'ctl' (charge-free dianion), 'bare' (charge-free monoanion), 'gctl'"""
    lmap = dict((l["key"], l) for l in LEVELS)
    if level in ("CAMSV", "CAMTZ"):
        src = lmap[level]["di"] if kind == "cond" else "%s_R/job.inp" % lmap[level]["ctl" if kind == "ctl" else "bare"]
        return open(os.path.join(s6, src)).read()
    tz = open(os.path.join(s6, lmap["CAMTZ"]["di"] if kind == "cond" else "d44/CAMTZ_C_ctl_R/job.inp")).read()
    if level == "CAMQZ":
        if tz.count("def2-TZVPD") != 1:
            stop("the committed CAMTZ input does not name its basis set exactly once.")
        return tz.replace("def2-TZVPD", "def2-QZVPPD")
    # CAMTZG: the committed geometry written out, and one ghost atom after it
    ghost = [g for g in sel["gctl"] if g["name"] == gname][0]["ghost"] if kind == "gctl" else cond["ghost"]
    ls = tz.split("\n")
    if ls[-2] != "* xyzfile %s substrate.xyz" % DI or ls[-1] != "":
        stop("the committed CAMTZ input does not end as this script expects.")
    geo = open(os.path.join(s6, GEOM[state])).read().split("\n")
    n = int(geo[0].split()[0])
    return "\n".join(ls[:-2] + ["* xyz %s" % DI] + geo[2:2 + n] + ["H: %.6f %.6f %.6f" % ghost, "*", ""])


def runs_of(sel):
    """every run in the order it is made within its job: dict(run, job, level, kind, cond, state, vpot, ref)"""
    out = []

    def add(run, job, level, kind, cond=None, state="R", vpot=False, ref=None, gname=None):
        out.append(dict(run=run, job=job, level=level, kind=kind, cond=cond, state=state, vpot=vpot, ref=ref, gname=gname))

    cs = sel["conds"]
    by = lambda part: [x for x in cs if x["part"] == part]
    # S1: the four controls with their maps, then every CAMSV run
    for kind, tag in (("ctl", "C_ctl"), ("bare", "C_bare")):
        for st in STATES:
            add("S1_%s_%s" % (tag, st), "S1", "CAMSV", kind, state=st, vpot=True)
    for x in by("A") + by("B") + by("C"):
        for st in x["levels"]["CAMSV"]:
            add(rname("CAMSV", x["name"], st), "S1", "CAMSV", "cond", x["name"], st)
    # T1: the four controls with their maps, then the ladder
    for kind, tag in (("ctl", "C_ctl"), ("bare", "C_bare")):
        for st in STATES:
            add("T1_%s_%s" % (tag, st), "T1", "CAMTZ", kind, state=st, vpot=True)
    for x in by("A"):
        for st in x["levels"]["CAMTZ"]:
            add(rname("CAMTZ", x["name"], st), "T1", "CAMTZ", "cond", x["name"], st)
    # T2: its control, the smaller size of part B first, then the rungs from the farthest in, then the larger size
    add("T2_C_ctl_R", "T2", "CAMTZ", "ctl")
    b = by("B")
    for x in sorted(b, key=lambda x: x["q"])[:len(b) // 2] + list(reversed(by("C"))) + sorted(b, key=lambda x: x["q"])[len(b) // 2:]:
        for st in x["levels"]["CAMTZ"]:
            add(rname("CAMTZ", x["name"], st), "T2", "CAMTZ", "cond", x["name"], st)
    # G1: a plain control, the ghost controls, then the ghost conditions, the smallest charge first
    add("G1_C_ctl_R", "G1", "CAMTZ", "ctl")
    for g in sel["gctl"]:
        for st in g["states"]:
            add("CAMTZG_%s_di_%s" % (g["name"], st), "G1", "CAMTZG", "gctl", state=st, gname=g["name"])
    for x in sorted([x for x in cs if "CAMTZG" in x["levels"]], key=lambda x: (x["q"], x["name"])):
        for st in x["levels"]["CAMTZG"]:
            add(rname("CAMTZG", x["name"], st), "G1", "CAMTZG", "cond", x["name"], st, ref="CAMTZG_%s_di_%s" % (x["gctl"], st))
    # Q1: state by state, the control first
    for st in STATES:
        add("Q1_C_ctl_%s" % st, "Q1", "CAMQZ", "ctl", state=st)
        for x in by("Q"):
            add(rname("CAMQZ", x["name"], st), "Q1", "CAMQZ", "cond", x["name"], st, ref="Q1_C_ctl_%s" % st)
    if len(set(r["run"] for r in out)) != len(out):
        stop("a run is listed twice.")
    return out


def points_text(sel):
    return "%d\n" % len(sel["points"]) + "".join(pt_line(p) + "\n" for p in sel["points"])


def conditions_text(sel):
    L = ["condition\tpart\twhere\tcharges\tsite\tq_e\tx\ty\tz\ta_kcal_per_e\tfirst_order_kcal\tnearest_atom\tnearest_state\tnearest_A\tpoint_row"
         "\tlevels_and_states\tghost_x\tghost_y\tghost_z\tghost_control\tpc_sha256"]
    f = lambda v, fmt: "" if v is None else fmt % v
    for x in sel["conds"]:
        g = x["ghost"]
        L.append("\t".join([x["name"], x["part"], x["where"], str(len(x["charges"])), f(x["site"], "%d"), f(x["q"], "%+.4f"),
                            f(x["xyz"] and x["xyz"][0], "%.6f"), f(x["xyz"] and x["xyz"][1], "%.6f"), f(x["xyz"] and x["xyz"][2], "%.6f"),
                            f(x["a"], "%+.6f"), f(x["first"], "%+.6f"), x["atom"] or "", x["kst"] or "", f(x["d"], "%.4f"), f(x["pt"], "%d"),
                            ";".join("%s:%s" % (k, ",".join(v)) for k, v in sorted(x["levels"].items())),
                            f(g and g[0], "%.6f"), f(g and g[1], "%.6f"), f(g and g[2], "%.6f"), x.get("gctl") or "",
                            hashlib.sha256(x["pc"].encode()).hexdigest()]))
    return "\n".join(L) + "\n"


def runs_text(sel):
    L = ["run\tjob\tlevel\tkind\tcondition\tstate\tmap\treference_control"]
    lmap = dict((l["key"], l) for l in LEVELS)
    for r in runs_of(sel):
        ref = r["ref"] or ("" if r["kind"] != "cond" else "%s_%s (committed)" % (lmap[r["level"]]["ctl"], r["state"]))
        L.append("\t".join([r["run"], r["job"], r["level"], r["kind"], r["cond"] or "", r["state"], "yes" if r["vpot"] else "no", ref]))
    return "\n".join(L) + "\n"


def levels_text(sel):
    L = ["key\tlevel\trole\tmethod_line\tinputs\tjobs"]
    how = {"CAMSV": "committed inputs, byte for byte", "CAMTZ": "committed inputs, byte for byte",
           "CAMTZG": "the committed CAMTZ input with the committed geometry written out and one ghost atom H: after it",
           "CAMQZ": "the committed CAMTZ input with def2-TZVPD changed to def2-QZVPPD"}
    for l in LEVELS:
        L.append("\t".join([l["key"], l["name"], l["role"], LINE[l["key"]], how[l["key"]], ",".join(j[0] for j in JOBS if j[1] == l["key"] or (l["key"] == "CAMTZ" and j[0] == "G1"))]))
    return "\n".join(L) + "\n"


def sources(sel):
    out = [GEOM["R"], GEOM["TS"], SITES, AMAT, FLOORS, D33A_PTS, D33A_VPOT["R"], D33A_VPOT["TS"]] + OTHER
    for lv in LEVELS[:2]:
        out.append(lv["di"])
        for c in ("ctl", "bare"):
            for st in STATES:
                out += ["%s_%s/job.inp" % (lv[c], st), "%s_%s/job.out" % (lv[c], st)]
    for d in DESIGNS_Q:
        out += ["d43/%s_di_%s/%s" % (d, st, f) for st in STATES for f in ("pc.pc", "job.out")]
    for d in D44_DESIGNS:
        out += ["d44/CAMTZ_%s_di_%s/%s" % (d, st, f) for st in STATES for f in ("pc.pc", "job.out")]
    for x in sel["conds"]:
        if x["part"] == "A" and abs(x["q"] - 1.0) < 1e-9:
            out += ["d43/%s_di_%s/%s" % (x["name"], st, f) for st in STATES for f in ("pc.pc", "job.out")]
    seen, uniq = set(), []
    for f in out:
        if f not in seen:
            seen.add(f)
            uniq.append(f)
    return uniq


def pbs_text(job, sel, runs, hpc):
    jk, lk, pbs = job
    mine = [r for r in runs if r["job"] == jk]
    L = ["#!/bin/bash", "#PBS -N %s" % pbs,
         "#PBS -l select=1:ncpus=8:mem=%s" % MEM, "#PBS -l walltime=%s" % WALL,
         "#PBS -m ae", "#PBS -M 18660916@sun.ac.za", "#PBS -j oe",
         "#PBS -o %s/run_%s.pbs.out" % (hpc, jk),
         "# Stage 6 DIAGNOSTIC 47, a preliminary: when does a bare positive charge take electrons from the substrate?",
         "# This job is one of %d. It makes %d runs in the order below; a run that fails does not stop the ones after it." % (len(JOBS), len(mine)),
         "# No -q: the allocation is automatic and walltime %s routes to day. Restartable:" % WALL,
         "#   submitted again, it skips every run that ended normally and every map that is written.",
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
         '  else ( cd "$d" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null ); '
         'echo "$1 orca_exit=$? normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$d/job.out") $(date)"; fi',
         "}",
         "vpot() {",
         '  d="$W/$1"; pts="$W/points.xyz"',
         '  if [ -s "$d/vpot.out" ]; then echo "$1 vpot done before"; return; fi',
         '  if ! grep -q "ORCA TERMINATED NORMALLY" "$d/job.out" 2>/dev/null; then echo "$1 vpot skipped, run not normal"; return; fi',
         '  ( cd "$d" && "$ORCA/orca_vpot" job.gbw job.scfp "$pts" vpot.out job > vpot.log 2>&1 )',
         '  echo "$1 vpot_exit=$? lines=$(wc -l < "$d/vpot.out" 2>/dev/null || echo 0) expected=%d $(date)"' % (len(sel["points"]) + 1),
         "}"]
    for r in mine:
        L.append("run %s" % r["run"])
        if r["vpot"]:
            L.append("vpot %s" % r["run"])
    L += ['echo "end=$(date)"',
          'n=0; for r in %s; do grep -q \'ORCA TERMINATED NORMALLY\' "$W/$r/job.out" 2>/dev/null && n=$((n+1)); done' % " ".join(r["run"] for r in mine),
          'echo "finished: $n of %d"' % len(mine)]
    return "\n".join(L) + "\n"


def nsum(sel, runs):
    n = 0
    for r in runs:
        n += 1                                                    # job.inp
        n += 0 if r["level"] == "CAMTZG" else 1                   # substrate.xyz
        n += 1 if r["kind"] == "cond" else 0                      # pc.pc
    return n + 5 + len(JOBS)                                      # CONDITIONS, RUNS, LEVELS, points.xyz, SOURCES; the job scripts


def main():
    if len(sys.argv) != 4:
        stop("usage: cs6_d47_prepare.py <results-repo> <out-dir> <hpc-dir>")
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
            f.write(input_text(s6, sel, r["level"], r["kind"], x, r["state"], r["gname"]))
        if r["level"] != "CAMTZG":
            shutil.copyfile(os.path.join(s6, GEOM[r["state"]]), os.path.join(d, "substrate.xyz"))
        if x is not None:
            with open(os.path.join(d, "pc.pc"), "w") as f:
                f.write(x["pc"])
        for fn in sorted(os.listdir(d)):
            sums.append((sha(os.path.join(d, fn)), "%s/%s" % (r["run"], fn)))
    for fn, text in (("CONDITIONS.tsv", conditions_text(sel)), ("RUNS.tsv", runs_text(sel)), ("LEVELS.tsv", levels_text(sel)), ("points.xyz", points_text(sel))):
        with open(os.path.join(out, fn), "w") as f:
            f.write(text)
    with open(os.path.join(out, "SOURCES.sha256"), "w") as f:
        for s in sources(sel):
            f.write("%s  %s\n" % (sha(os.path.join(s6, s)), s))
    for fn in ("CONDITIONS.tsv", "RUNS.tsv", "LEVELS.tsv", "points.xyz", "SOURCES.sha256"):
        sums.append((sha(os.path.join(out, fn)), fn))
    for job in JOBS:
        p = os.path.join(out, "run_%s.pbs" % job[0])
        with open(p, "w") as f:
            f.write(pbs_text(job, sel, runs, hpc))
        os.chmod(p, 0o755)
        sums.append((sha(p), "run_%s.pbs" % job[0]))
    if len(sums) != nsum(sel, runs):
        stop("%d files were written, not %d." % (len(sums), nsum(sel, runs)))
    # the first NSITE lines of the points file are d33a's, byte for byte
    if open(os.path.join(out, "points.xyz")).read().split("\n")[1:1 + NSITE] != open(os.path.join(s6, D33A_PTS)).read().split("\n")[1:1 + NSITE]:
        stop("the points file does not begin with the %d lines of %s." % (NSITE, D33A_PTS))
    with open(os.path.join(out, "INPUT_CHECKSUMS.sha256"), "w") as f:
        for h, n in sums:
            f.write("%s  %s\n" % (h, n))
    for x in sel["conds"]:
        print("   %-22s part %s  %-40s  %s" % (x["name"], x["part"], x["where"], ";".join("%s:%s" % (k, ",".join(v)) for k, v in sorted(x["levels"].items()))))
    for jk, lk, pbs in JOBS:
        print("   job %s  %-7s %2d runs" % (jk, lk, sum(1 for r in runs if r["job"] == jk)))
    print("   %d conditions, %d run directories, %d files checksummed, %d jobs" % (len(sel["conds"]), len(runs), len(sums), len(JOBS)))


if __name__ == "__main__":
    try:
        main()
    except Refused as e:
        sys.stderr.write("REFUSED  %s\n" % e)
        raise SystemExit(2)
