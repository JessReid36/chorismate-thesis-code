#!/usr/bin/env python3
"""
s18_provenance_check.py - provenance and consistency audit of the ensemble files that
Phase 2.2 rests on.

Reads ONLY committed files in chorismate-thesis-results. Changes nothing. Prints one
line per check, each marked

    PASS   the check holds
    FAIL   the check does not hold, for something Phase 2.2 depends on
    GAP    an artefact that should exist for provenance is not committed
    INFO   a fact worth recording that is not itself a pass/fail test

and exits non-zero if any check FAILS.

USAGE
    python3 s18_provenance_check.py <results_repo> [<code_repo>]

The code repo is optional. With it, section 7 imports s16's own Kabsch fit and CORE
list, so the grid check uses exactly the alignment that built the grid, and section 5
checks the atom indices hard-coded in s11, s13 and s16 against the chemistry.

SECTIONS
    1  frame sets: which frames exist, and which manifest defines them
    2  inventory: every per-frame artefact, by stage, for post-cut and pilot frames
    3  inputs: level of theory, SCF settings, charge, QM region, restraints
    4  outputs: the 20_invacuo single points
    5  atom order, indexing and connectivity
    6  geometry identity across every committed copy of the same structure
    7  grid_v2 against all 30 aligned frames, and A_v2's frame order
    8  s17 re-derived from raw energies
    9  Phase 1 barrier provenance: ensemble_barriers.tsv against path_summary.txt

Python 3.6 compatible. Standard library except section 7, which needs numpy.
"""
import math
import re
import sys
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path

H2KCAL = 627.5094740631          # the factor s17 uses
CANON_NAMES = ("C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 "
               "C10 O5 O6").split()
QM_FIRST, QM_LAST = 6208, 6231   # 1-based rows / PDB serials of CHA 383
LOT = {"B3LYP", "D3BJ", "DEF2-SVP", "DEF2/J", "RIJCOSX"}

counts = Counter()


def say(tag, msg):
    counts[tag] += 1
    print("  {:<5} {}".format(tag, msg))


def head(title):
    print("\n" + "=" * 78 + "\n" + title + "\n" + "=" * 78)


def read(p):
    return Path(p).read_text(errors="replace")


def read_xyz(p):
    L = read(p).splitlines()
    n = int(L[0].split()[0])
    els, xyz = [], []
    for l in L[2:2 + n]:
        f = l.split()
        els.append(f[0])
        xyz.append(tuple(float(v) for v in f[1:4]))
    return els, xyz


def maxdiff(a, b):
    return max(abs(x - y) for p, q in zip(a, b) for x, y in zip(p, q))


def dist(p, q):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(p, q)))


def out_coords(text):
    m = re.search(r"CARTESIAN COORDINATES \(ANGSTROEM\)\n-+\n(.*?)\n\s*\n", text, re.S)
    if not m:
        return None, None
    els, xyz = [], []
    for l in m.group(1).splitlines():
        f = l.split()
        if len(f) >= 4:
            els.append(f[0])
            xyz.append(tuple(float(v) for v in f[1:4]))
    return els, xyz


def last_energy(text):
    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    return (float(m[-1]) if m else None), len(m)


def read_tsv(p):
    hdr, rows = None, []
    for l in read(p).splitlines():
        if not l.strip():
            continue
        if l.startswith("#"):
            # a commented header row, e.g. "# idx<TAB>frame<TAB>...", is still the header
            if hdr is None and "\t" in l:
                hdr = l.lstrip("#").strip().split("\t")
            continue
        f = l.split("\t")
        if hdr is None:
            hdr = f
            continue
        rows.append(dict(zip(hdr, f)))
    return hdr, rows


# ------------------------------------------------------------------ connectivity
RCOV = {"H": 0.31, "C": 0.76, "O": 0.66}


def bonds(els, xyz, scale=1.25):
    b = set()
    for i in range(len(els)):
        for j in range(i + 1, len(els)):
            if dist(xyz[i], xyz[j]) < scale * (RCOV[els[i]] + RCOV[els[j]]):
                b.add((i, j))
    return frozenset(b)


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    R = Path(sys.argv[1])
    C = Path(sys.argv[2]) if len(sys.argv) > 2 else None
    Q = R / "05_qmmm"
    INV, ENS, HAR = Q / "20_invacuo", Q / "19_ensemble_barriers", Q / "19_ensemble_harvest"
    SEL = Q / "12_frame_selection"
    for p in (INV, ENS, Q / "13_bridge" / "complex_solvated.ORCAFF.prms",
              R / "phase2.2" / "grid_v2.tsv", R / "phase2.2" / "A_v2.tsv"):
        if not p.exists():
            sys.exit("FAIL: required input missing: {}".format(p))
    STATES = OrderedDict([("R", "reactant"), ("TS", "transition_state"), ("P", "product")])

    # ============================================================== 1 frame sets
    head("1. FRAME SETS")
    _, tsv = read_tsv(ENS / "ensemble_barriers.tsv")
    tsv_frames = [r["frame"] for r in tsv]
    say("INFO", "ensemble_barriers.tsv: {} frames".format(len(tsv_frames)))
    gl = [l for l in read(R / "phase2.2" / "grid_v2.tsv").splitlines()
          if l.startswith("#   frames =")]
    post = gl[0].split("=", 1)[1].strip().split(",")
    pilot = [f for f in tsv_frames if f not in post]
    say("PASS" if all(f in tsv_frames for f in post) else "FAIL",
        "all {} grid_v2 frames are in ensemble_barriers.tsv".format(len(post)))
    say("INFO", "post-cut (grid_v2): {}   pilot: {} ({})".format(
        len(post), len(pilot), " ".join(sorted(pilot))))
    _, ext = read_tsv(SEL / "selection_manifest_ext_only.tsv")
    _, ren = read_tsv(SEL / "selection_manifest_ext_renumbered.tsv")
    ok = len(ext) == len(ren) and all(
        int(b["frame"]) == int(a["frame"]) + 20000 and
        abs(float(b["prod_ps"]) - float(a["prod_ps"]) - 20000.0) < 1e-6 and
        all(a[k] == b[k] for k in a if k not in ("frame", "prod_ps"))
        for a, b in zip(ext, ren))
    say("PASS" if ok else "FAIL", "ext_renumbered = ext_only with frame and prod_ps + 20000 "
        "(the extension trajectory starts at 20 ns), every other column identical "
        "({} rows)".format(len(ren)))
    ren_frames = ["{:05d}".format(int(r["frame"])) for r in ren]
    say("PASS" if sorted(ren_frames) == sorted(post) else "FAIL",
        "the 30 grid_v2 frames are exactly selection_manifest_ext_renumbered")
    sel = []
    for name in ("selection_manifest.tsv", "selection_manifest_nac_extra.tsv"):
        _, rows = read_tsv(SEL / name)
        sel += ["{:05d}".format(int(r["frame"])) for r in rows]
    say("PASS" if sorted(set(pilot)) == sorted(set(sel) & set(pilot)) else "FAIL",
        "every pilot frame is in selection_manifest(.tsv | _nac_extra.tsv)")
    unused = sorted(set(sel) - set(tsv_frames))
    if unused:
        say("INFO", "selected but never run through the ensemble: {}".format(" ".join(unused)))
    say("INFO", "selection_manifest_eq30.tsv (frames 11207-19999) is not used by any "
        "committed ensemble frame")

    # ============================================================== 2 inventory
    head("2. INVENTORY: per-frame artefacts by stage")
    art = OrderedDict()
    art["12 frame rst7"] = lambda f: SEL / "frames" / "frame_{}_CHA2.rst7".format(f)
    art["12 frame qmatoms"] = lambda f: SEL / "frames" / "frame_{}_CHA2.qmatoms".format(f)
    art["19h reactant_opt.inp"] = lambda f: HAR / "frame_{}".format(f) / "reactant_opt.inp"
    art["19h reactant_qm.xyz"] = lambda f: HAR / "frame_{}".format(f) / "reactant_qm.xyz"
    art["19 neb.inp"] = lambda f: ENS / "frame_{}".format(f) / "neb.inp"
    art["19 path_summary.txt"] = lambda f: ENS / "frame_{}".format(f) / "path_summary.txt"
    for s, long in STATES.items():
        art["19 {}_qm.xyz".format(long)] = (lambda L: lambda f: ENS / "frame_{}".format(f)
                                            / "{}_qm.xyz".format(L))(long)
        art["19 invacuo_{}.inp".format(long)] = (lambda L: lambda f: ENS / "frame_{}".format(f)
                                                 / "invacuo_{}.inp".format(L))(long)
    for s in STATES:
        for ext_ in ("xyz",):
            art["20 {}_{}.xyz".format("<f>", s)] = (lambda S: lambda f: INV / "{}_{}.xyz".format(f, S))(s)
        for ext_ in ("inp", "out", "property.txt"):
            art["20 sp_<f>_{}.{}".format(s, ext_)] = (lambda S, E: lambda f: INV / "sp_{}_{}.{}".format(f, S, E))(s, ext_)
    needed = {"20 <f>_R.xyz", "20 <f>_TS.xyz", "20 sp_<f>_R.out", "20 sp_<f>_TS.out"}
    print("  {:<26}{:>10}{:>9}   post-cut frames missing".format("artefact", "post-cut", "pilot"))
    gaps = []
    for name, fn in art.items():
        mp = [f for f in post if not fn(f).exists()]
        mq = [f for f in pilot if not fn(f).exists()]
        miss = " ".join(mp) if len(mp) <= 8 else "{} frames".format(len(mp))
        print("  {:<26}{:>7}/{:<2}{:>6}/{:<2}   {}".format(
            name, len(post) - len(mp), len(post), len(pilot) - len(mq), len(pilot), miss))
        if mp:
            gaps.append((name, len(mp)))
            if name in needed:
                say("FAIL", "{}: missing for post-cut frames {}".format(name, " ".join(mp)))
    for name, n in gaps:
        if name not in needed:
            say("GAP", "{} not committed for {} of {} post-cut frames".format(name, n, len(post)))
    say("GAP", "no per-frame scan inputs or outputs are committed for any ensemble frame. "
        "The scan restraint (run_scan.sh: harmonic colvar restraints, SPRING in kJ/mol/A^2; "
        "2500 now, 400 for the pilot frames per the phase2.2 README) cannot be verified per "
        "frame from the repository")
    say("INFO", "frame 08170 (pilot) has no TS geometry anywhere: 19 and 20 both lack it")

    # ============================================================== 3 inputs
    head("3. INPUTS: level of theory, settings, charge, QM region, restraints")

    def simple_tokens(text):
        toks = set()
        for l in text.splitlines():
            if l.strip().startswith("!"):
                toks |= {t.upper() for t in l.strip()[1:].split()}
        return toks

    def coord_line(text):
        m = re.search(r"^\*\s*(xyzfile|pdbfile)\s+(-?\d+)\s+(\d+)\s+(\S+)", text, re.M | re.I)
        return m.groups() if m else None

    groups = OrderedDict()
    groups["20 sp_*.inp"] = sorted(INV.glob("sp_*.inp"))
    groups["19 invacuo_*.inp"] = sorted(ENS.glob("frame_*/invacuo_*.inp"))
    groups["19 neb.inp"] = sorted(ENS.glob("frame_*/neb.inp"))
    groups["19h reactant_opt.inp"] = sorted(HAR.glob("frame_*/reactant_opt.inp"))
    for g, files in groups.items():
        texts = {p: read(p) for p in files}
        bad = [p for p, t in texts.items() if not LOT <= simple_tokens(t)]
        extras = Counter(" ".join(sorted(simple_tokens(t) - LOT)) for t in texts.values())
        say("PASS" if not bad else "FAIL",
            "{} ({} files): B3LYP D3BJ def2-SVP def2/J RIJCOSX in every file".format(g, len(files)))
        say("INFO", "{} other simple-input keywords: {}".format(
            g, "; ".join("[{}] x{}".format(k, v) for k, v in extras.items())))
        cl = Counter(coord_line(t)[:3] if coord_line(t) else None for t in texts.values())
        say("PASS" if set(k[1:] for k in cl if k) == {("-2", "1")} and None not in cl else "FAIL",
            "{}: charge -2, multiplicity 1 in every file {}".format(g, dict(cl)))
        blocks = Counter(tuple(sorted(set(re.findall(r"^%(\w+)", t, re.M | re.I))))
                         for t in texts.values())
        say("INFO", "{} % blocks: {}".format(g, "; ".join("{} x{}".format(" ".join(k), v)
                                                          for k, v in blocks.items())))
        rst = [p for p, t in texts.items()
               if re.search(r"constrain|restrain|colvar|\bscan\b|fixed", t, re.I)]
        say("INFO", "{}: {} file(s) contain constraint/restraint/scan keywords".format(g, len(rst)))

    # exact-text variants after masking the frame-specific tokens
    def physics(t):
        """the input with run-control text removed: %pal lines and the KeepDens keyword"""
        L = [l for l in t.splitlines() if not re.match(r"^\s*%pal\b", l, re.I)]
        L = [re.sub(r"\s+KeepDens\b", "", l, flags=re.I) for l in L]
        return "\n".join(l.rstrip() for l in L if l.strip())

    def mask_sp(p, t):
        f, s = re.match(r"sp_(\d+)_(\w+)\.inp", p.name).groups()
        own = "{}_{}.xyz".format(f, s)
        return t.replace(own, "<f>_<s>.xyz"), t.count(own)

    var, phys = Counter(), Counter()
    wrong = []
    for p in groups["20 sp_*.inp"]:
        m, n = mask_sp(p, read(p))
        var[m] += 1
        phys[physics(m)] += 1
        if n != 1:
            wrong.append(p.name)
    say("PASS" if len(phys) == 1 else "FAIL",
        "20 sp_*.inp: {} variant(s) once run-control text (%pal, KeepDens) is set aside".format(len(phys)))
    say("INFO", "20 sp_*.inp: {} literal variants, differing only in %pal (parallel ranks) "
        "and KeepDens (file retention); neither changes the result".format(len(var)))
    say("PASS" if not wrong else "FAIL",
        "20 sp_*.inp: every input reads its own <frame>_<state>.xyz {}".format(wrong or ""))
    same, diff_ = 0, []
    for p in groups["19 invacuo_*.inp"]:
        f = p.parent.name.replace("frame_", "")
        s = {v: k for k, v in STATES.items()}[p.stem.replace("invacuo_", "")]
        q = INV / "sp_{}_{}.inp".format(f, s)
        if q.exists() and read(q) == read(p):
            same += 1
        else:
            diff_.append(p.parent.name + "/" + p.name)
    say("PASS" if not diff_ else "FAIL",
        "19 invacuo_*.inp byte-identical to 20 sp_*.inp: {}/{} {}".format(
            same, same + len(diff_), diff_[:5] or ""))

    def split_active(t):
        m = re.search(r"ActiveAtoms\s*\{([^}]*)\}", t)
        act = set(int(x) for x in m.group(1).split()) if m else set()
        return re.sub(r"ActiveAtoms\s*\{[^}]*\}", "ActiveAtoms {<list>}", t), act

    for g in ("19 neb.inp", "19h reactant_opt.inp"):
        rest, acts = Counter(), Counter()
        qm_ok = True
        for p in groups[g]:
            t = read(p).replace(p.parent.name.replace("frame_", ""), "<f>")
            r_, a_ = split_active(t)
            rest[r_] += 1
            acts[frozenset(a_)] += 1
            qm = re.search(r"QMAtoms\s*\{(\d+):(\d+)\}", t)
            if not qm or (int(qm.group(1)), int(qm.group(2))) != (QM_FIRST - 1, QM_LAST - 1):
                qm_ok = False
            elif not set(range(QM_FIRST - 1, QM_LAST)) <= a_:
                qm_ok = False
        if g == "19 neb.inp":
            kind = defaultdict(list)
            for p in groups[g]:
                k = "NEB-TS" if re.search(r"\bNEB-TS\b", read(p), re.I) else (
                    "NEB-CI" if re.search(r"\bNEB-CI\b", read(p), re.I) else "other")
                kind[k].append(p.parent.name.replace("frame_", ""))
            rest2 = Counter(re.sub(r"\bNEB-(TS|CI)\b", "NEB-<x>", r_) for r_ in rest.elements())
            say("PASS" if len(rest2) == 1 else "FAIL",
                "{}: {} variant(s) outside ActiveAtoms once NEB-TS/NEB-CI is set aside".format(g, len(rest2)))
            for k, fs in kind.items():
                say("INFO", "{} {} x{}: {}".format(g, k, len(fs), " ".join(fs)))
        else:
            say("PASS" if len(rest) == 1 else "FAIL",
                "{}: {} text variant(s) outside the ActiveAtoms list, frame number masked".format(g, len(rest)))
        say("INFO", "{}: {} distinct ActiveAtoms set(s), sizes {}".format(
            g, len(acts), sorted({len(a) for a in acts})))
        say("PASS" if qm_ok else "FAIL",
            "{}: QMAtoms {{{}:{}}} (0-based) in every file, inside ActiveAtoms".format(
                g, QM_FIRST - 1, QM_LAST - 1))
    nt = simple_tokens(read(groups["19 neb.inp"][0])) - LOT
    ot = simple_tokens(read(groups["19h reactant_opt.inp"][0])) - LOT if groups["19h reactant_opt.inp"] else set()
    say("INFO", "same functional, dispersion, basis and RI as the single points at every "
        "stage; SCF settings differ by stage (neb: {}; reactant opt: {})".format(
            " ".join(sorted(nt)), " ".join(sorted(ot))))

    # ============================================================== 4 outputs
    head("4. OUTPUTS: 20_invacuo single points")
    outs = sorted(INV.glob("sp_*.out"))
    prop = defaultdict(Counter)
    problems = defaultdict(list)
    E = {}
    for p in outs:
        t = read(p)
        key = p.stem.replace("sp_", "")
        v = re.search(r"Program Version (\S+)", t)
        prop["version"][v.group(1) if v else None] += 1
        for lab, pat in (("charge", r"Total Charge\s+Charge\s+\.+\s+(-?\d+)"),
                         ("mult", r"Multiplicity\s+Mult\s+\.+\s+(\d+)"),
                         ("NEL", r"Number of Electrons\s+NEL\s+\.+\s+(\d+)"),
                         ("basis dim", r"Basis Dimension\s+Dim\s+\.+\s+(\d+)")):
            m = re.search(pat, t)
            prop[lab][m.group(1) if m else None] += 1
        if "TERMINATED NORMALLY" not in t:
            problems["not terminated normally"].append(key)
        if not re.search(r"SCF CONVERGED AFTER\s+\d+\s+CYCLES", t) or "SCF NOT CONVERGED" in t:
            problems["SCF not converged"].append(key)
        e, n = last_energy(t)
        if n != 1:
            problems["FINAL SINGLE POINT ENERGY count != 1"].append(key)
        E[key] = e
        echo = [re.sub(r"^\|\s*\d+>\s?", "", l).rstrip()
                for l in t.splitlines() if re.match(r"^\|\s*\d+>", l)]
        echo = [l for l in echo if l.strip() and "END OF INPUT" not in l]
        inp = [l.rstrip() for l in read(p.with_suffix(".inp")).splitlines() if l.strip()]
        if echo != inp:
            problems["input echo differs from the committed .inp"].append(key)
        w = re.search(r"WARNINGS\n.*?carefully!\n=+\n(.*?)\n=+\n", t, re.S)
        if w and w.group(1).strip():
            problems["ORCA WARNINGS block not empty"].append(key)
    expect = {"version": "6.0.1", "charge": "-2", "mult": "1", "NEL": "118"}
    for lab, c in prop.items():
        ok = len(c) == 1 and (lab not in expect or list(c)[0] == expect[lab])
        say("PASS" if ok else "FAIL", "{:<10} {} across {} outputs".format(lab, dict(c), len(outs)))
    for lab in ("not terminated normally", "SCF not converged",
                "FINAL SINGLE POINT ENERGY count != 1",
                "input echo differs from the committed .inp", "ORCA WARNINGS block not empty"):
        say("PASS" if not problems[lab] else "FAIL",
            "{}: {} {}".format(lab, len(problems[lab]), problems[lab][:6] or ""))
    worst, nmatch, missing = 0.0, 0, []
    for r in tsv:
        for s, col in (("R", "E_vac_R"), ("TS", "E_vac_TS"), ("P", "E_vac_P")):
            k = "{}_{}".format(r["frame"], s)
            if not r.get(col):
                continue
            if E.get(k) is None:
                missing.append(k)
                continue
            worst = max(worst, abs(E[k] - float(r[col])))
            nmatch += 1
    say("PASS" if worst <= 5.1e-7 and not missing else "FAIL",
        "E_vac_R/TS/P in ensemble_barriers.tsv = last FINAL SINGLE POINT ENERGY: {} values, "
        "max |diff| {:.1e} Eh (tsv prints 6 dp) {}".format(nmatch, worst, missing or ""))
    bw = 0.0
    for r in tsv:
        k1, k2 = r["frame"] + "_R", r["frame"] + "_TS"
        if r.get("barrier_vac") and E.get(k1) is not None and E.get(k2) is not None:
            bw = max(bw, abs((E[k2] - E[k1]) * H2KCAL - float(r["barrier_vac"])))
    say("PASS" if bw < 1e-5 else "FAIL",
        "barrier_vac = (E_TS - E_R) x {}: max |diff| {:.1e} kcal/mol".format(H2KCAL, bw))

    # ============================================================== 5 indexing
    head("5. ATOM ORDER, INDEXING AND CONNECTIVITY")
    L = read(Q / "13_bridge" / "complex_solvated.ORCAFF.prms").splitlines()
    natom = int(L[3].split()[0])
    rows = [l.split() for l in L[4:4 + natom]]
    sub = [r for r in rows if QM_FIRST <= int(r[0]) <= QM_LAST]
    canon = [r[1] for r in sub]
    qsum = sum(float(r[2]) for r in sub)
    say("PASS" if len(sub) == 24 and abs(qsum + 2) < 1e-6 else "FAIL",
        "ORCAFF.prms rows {}-{}: 24 atoms, charge {:+.6f}, order {}".format(
            QM_FIRST, QM_LAST, qsum, "".join(canon)))
    pdb = Q / "18e_reactant_reduced" / "reactant_reduced.pdb"
    if pdb.exists():
        pa = [(int(l[6:11]), l[12:16].strip(), l[76:78].strip() or l[12:16].strip()[0])
              for l in read(pdb).splitlines()
              if l.startswith(("ATOM", "HETATM")) and l[17:20] == "CHA" and int(l[22:26]) == 383]
        ok = ([a[0] for a in pa] == list(range(QM_FIRST, QM_LAST + 1)) and
              [a[1] for a in pa] == CANON_NAMES and [a[2] for a in pa] == canon)
        say("PASS" if ok else "FAIL",
            "PDB CHA 383: serials {}-{}, names {} ... {}, elements match the prms".format(
                QM_FIRST, QM_LAST, " ".join(CANON_NAMES[:8]), CANON_NAMES[-1]))
    qf = sorted((SEL / "frames").glob("*.qmatoms"))
    bad = [p.name for p in qf if [int(x) for x in read(p).split()] != list(range(QM_FIRST, QM_LAST + 1))]
    say("PASS" if not bad else "FAIL", "{} .qmatoms files list {}-{} (1-based) {}".format(
        len(qf), QM_FIRST, QM_LAST, bad or ""))
    rs = sorted((SEL / "frames").glob("*.rst7"))
    bad = [p.name for p in rs if int(read(p).splitlines()[1].split()[0]) != natom]
    say("PASS" if not bad else "FAIL", "{} .rst7 files hold {} atoms, as the prms {}".format(
        len(rs), natom, bad or ""))

    xyzs = sorted(INV.glob("[0-9]*_*.xyz")) + sorted(ENS.glob("frame_*/*_qm.xyz")) + \
        sorted(HAR.glob("frame_*/reactant_qm.xyz"))
    bad = []
    geo = {}
    for p in xyzs:
        els, xyz = read_xyz(p)
        geo[p] = (els, xyz)
        if els != canon:
            bad.append(str(p.relative_to(R)))
    say("PASS" if not bad else "FAIL",
        "{} committed .xyz files: 24 atoms in the prms order {}".format(len(xyzs), bad[:5] or ""))
    tsc = Counter(read(p).splitlines()[1].strip() for p in xyzs
                  if p.name.endswith("_TS.xyz") or p.name == "transition_state_qm.xyz")
    say("PASS" if list(tsc) == ["Coordinates from ORCA-job neb_NEB-CI_converged.QMRegion"] else "FAIL",
        "every committed TS geometry is the converged climbing image {}, so NEB-TS and NEB-CI "
        "frames share one TS definition".format(dict(tsc)))

    ref = INV / "{}_R.xyz".format(post[0])
    rb = bonds(*geo[ref])
    say("INFO", "reference bond graph: {} bonds, from {}".format(len(rb), ref.name))
    nm = lambda b: "{}-{}".format(CANON_NAMES[b[0]], CANON_NAMES[b[1]])
    badR, change = [], defaultdict(Counter)
    for f in post + pilot:
        pR = INV / "{}_R.xyz".format(f)
        if not pR.exists():
            continue
        bR = bonds(*geo[pR])
        if bR != rb:
            badR.append("{}: +{} -{}".format(f, [nm(b) for b in bR - rb], [nm(b) for b in rb - bR]))
        for s in ("TS", "P"):
            pS = INV / "{}_{}.xyz".format(f, s)
            if pS.exists():
                bS = bonds(*geo[pS])
                key = "formed {} broken {}".format(sorted(nm(b) for b in bS - bR) or "none",
                                                   sorted(nm(b) for b in bR - bS) or "none")
                change[s][key] += 1
    say("PASS" if not badR else "FAIL",
        "every reactant (20_invacuo) has the reference bond graph {}".format(badR[:4] or ""))
    for s in ("TS", "P"):
        for k, v in change[s].items():
            say("INFO", "{} vs its own R, x{}: {}".format(s, v, k))
    say("PASS" if len(change["P"]) == 1 else "FAIL",
        "every product differs from its reactant by the same bond changes")

    # the indices the scripts hard-code
    els0, xyz0 = geo[ref]
    b0 = rb
    nb = lambda i: sorted(j for b in b0 for j in b if i in b and j != i)
    checks = [("index 7 is the ether O3, bonded to C2 (3) and C4 (8)",
               els0[7] == "O" and nb(7) == [3, 8]),
              ("indices 5, 6 are carboxylate O on C3 (4)", all(els0[i] == "O" and nb(i) == [4] for i in (5, 6))),
              ("indices 22, 23 are carboxylate O on C10 (21)", all(els0[i] == "O" and nb(i) == [21] for i in (22, 23))),
              ("index 19 is the hydroxyl O4 (bonded to C9 (17) and H8 (20))", els0[19] == "O" and nb(19) == [17, 20])]
    for lab, ok in checks:
        say("PASS" if ok else "FAIL", lab)
    if C:
        s16 = read(C / "phase2.2" / "s16_align_frames.py")
        core = [int(x) for x in re.search(r"^CORE = \[([^\]]*)\]", s16, re.M).group(1).split(",")]
        ring = {8, 10, 12, 13, 15, 17}
        say("PASS" if set(core) == ring | {3, 7} else "FAIL",
            "s16 CORE {} = ring carbons {} + O3 (7) + C2 (3); its comment calls it "
            "'ring carbons plus the ether oxygen'".format(core, sorted(ring)))
        for name in ("s11_probe_distance_scan.py", "s13_lj_charge_site_test.sh"):
            t = read(C / "phase2.2" / name)
            say("INFO", "{} refers to O3 as index 7: {}".format(
                name, bool(re.search(r"at\[7\]|O3_INDEX = 7|index 7", t))))

    # ============================================================== 6 geometry identity
    head("6. GEOMETRY IDENTITY ACROSS COMMITTED COPIES")
    w1, n1, b1 = 0.0, 0, []
    for f in post + pilot:
        for s, long in STATES.items():
            a, b = INV / "{}_{}.xyz".format(f, s), ENS / "frame_{}".format(f) / "{}_qm.xyz".format(long)
            if a.exists() and b.exists():
                d = maxdiff(geo[a][1], geo[b][1])
                w1, n1 = max(w1, d), n1 + 1
                if d > 1e-6:
                    b1.append("{}_{}".format(f, s))
    say("PASS" if not b1 else "FAIL",
        "20_invacuo <f>_<s>.xyz = 19 *_qm.xyz: {} pairs, max |dx| {:.1e} A {}".format(n1, w1, b1[:5] or ""))
    w2, n2, b2 = 0.0, 0, []
    for p in outs:
        key = p.stem.replace("sp_", "")
        g = INV / "{}.xyz".format(key)
        els, xyz = out_coords(read(p))
        if xyz is None or els != canon:
            b2.append(key)
            continue
        d = maxdiff(xyz, geo[g][1])
        w2, n2 = max(w2, d), n2 + 1
        if d > 1e-5:
            b2.append(key)
    say("PASS" if not b2 else "FAIL",
        "sp_*.out CARTESIAN COORDINATES = the .xyz it read: {} outputs, max |dx| {:.1e} A "
        "(ORCA prints 6 dp) {}".format(n2, w2, b2[:5] or ""))
    say("INFO", "so every density lives in the ORIGINAL QM/MM coordinates, not in "
        "s16-aligned space")
    w3, n3, b3 = 0.0, 0, []
    for p in sorted(HAR.glob("frame_*/reactant_qm.xyz")):
        q = ENS / p.parent.name / "reactant_qm.xyz"
        if q.exists():
            d = maxdiff(geo[p][1], geo[q][1])
            w3, n3 = max(w3, d), n3 + 1
            if d > 1e-6:
                b3.append(p.parent.name)
    say("PASS" if not b3 else "FAIL",
        "19_ensemble_harvest reactant = 19_ensemble_barriers reactant: {} pairs, max |dx| "
        "{:.1e} A {}".format(n3, w3, b3[:5] or ""))
    sx = R / "phase2.2" / "s13_lj_test2" / "sub.xyz"
    if sx.exists():
        d = maxdiff(read_xyz(sx)[1], geo[INV / "24883_R.xyz"][1])
        say("PASS" if d <= 1e-6 else "FAIL", "s13 sub.xyz = 24883_R.xyz, max |dx| {:.1e} A".format(d))

    # ============================================================== 7 grid
    head("7. GRID_v2 AGAINST ALL 30 ALIGNED FRAMES")
    ah = [l for l in read(R / "phase2.2" / "A_v2.tsv").splitlines() if l.startswith("# frames:")]
    afr = ah[0].split()[3:]
    say("PASS" if afr == post else "FAIL",
        "A_v2 frame columns are the grid_v2 frames in the same order")
    try:
        import numpy as np
    except ImportError:
        say("INFO", "numpy not available: section 7 skipped")
        np = None
    if np is not None:
        if C:
            import importlib.util
            spec = importlib.util.spec_from_file_location("s16", str(C / "phase2.2" / "s16_align_frames.py"))
            s16m = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(s16m)
            kabsch, core = s16m.kabsch, s16m.CORE
            say("INFO", "alignment: s16_align_frames.kabsch and CORE, imported from the code repo")
        else:
            core = [3, 7, 8, 10, 12, 13, 15, 17]

            def kabsch(P, Qm, idx):
                Pc, Qc = P[idx].mean(0), Qm[idx].mean(0)
                V, S, W = np.linalg.svd((P[idx] - Pc).T @ (Qm[idx] - Qc))
                d = np.sign(np.linalg.det(V @ W))
                Rm = V @ np.diag([1.0, 1.0, d]) @ W
                return Rm, Qc - Pc @ Rm
            say("INFO", "alignment: local copy of s16's Kabsch fit (no code repo given)")

        def geom(f, s):
            a = ENS / "frame_{}".format(f) / "{}_qm.xyz".format(STATES[s])
            b = INV / "{}_{}.xyz".format(f, s)
            return np.array(read_xyz(a if a.exists() else b)[1])
        G = np.array([[float(v) for v in l.split("\t")[1:4]]
                      for l in read(R / "phase2.2" / "grid_v2.tsv").splitlines()
                      if l.strip() and not l.startswith(("#", "idx"))])
        ref_xyz = geom(post[0], "R")
        atoms, perframe = [], []
        for f in post:
            Rm, t = kabsch(geom(f, "R"), ref_xyz, core)
            X = np.vstack([geom(f, "R") @ Rm + t, geom(f, "TS") @ Rm + t])
            atoms.append(X)
            perframe.append(np.linalg.norm(G[:, None, :] - X[None, :, :], axis=2).min())
        A = np.vstack(atoms)
        d = np.linalg.norm(G[:, None, :] - A[None, :, :], axis=2).min(1)
        say("PASS" if np.all(np.abs(d - 2.5) <= 1e-3) else "FAIL",
            "{} sites against the union of {} aligned frames (R+TS): min {:.4f} max {:.4f} A; "
            "{} at 2.500 +- 0.001".format(len(G), len(post), d.min(), d.max(),
                                          int((np.abs(d - 2.5) <= 1e-3).sum())))
        say("PASS" if min(perframe) >= 2.5 - 1e-4 else "FAIL",
            "no frame has any site closer than 2.500 A: smallest per-frame minimum {:.4f} A".format(
                min(perframe)))

    # ============================================================== 8 s17
    head("8. s17 RE-DERIVED FROM RAW ENERGIES")
    S17 = R / "phase2.2" / "s17_check"
    if not S17.exists():
        say("GAP", "phase2.2/s17_check not committed")
    else:
        tab = {}
        for l in read(S17 / "s17_check.pbs.out").splitlines():
            f = l.split()
            if len(f) == 5 and re.match(r"^\d{5}$", f[0]):
                tab[f[0]] = (float(f[1]), float(f[2]))
        Arows = [[float(x) for x in l.split("\t")] for l in read(R / "phase2.2" / "A_v2.tsv").splitlines()
                 if l.strip() and not l.startswith("#")]
        site = 117
        for fr, (lin_s, dft_s) in sorted(tab.items()):
            eR, eT = E.get(fr + "_R"), E.get(fr + "_TS")
            fR, _ = last_energy(read(S17 / "field_{}_{}_R.out".format(site, fr)))
            fT, _ = last_energy(read(S17 / "field_{}_{}_TS.out".format(site, fr)))
            dft = ((fT - eT) - (fR - eR)) * H2KCAL
            lin = Arows[site][afr.index(fr)]
            ok = abs(dft - dft_s) <= 1e-3 and abs(lin - lin_s) <= 1e-3
            say("PASS" if ok else "FAIL",
                "frame {}: direct DFT {:+.4f} (s17 {:+.3f}), linear {:+.4f} (s17 {:+.3f}), "
                "ratio {:.3f}".format(fr, dft, dft_s, lin, lin_s, dft / lin))
        echo_ok = True
        for p in sorted(S17.glob("field_*.out")):
            t = read(p)
            ec = [re.sub(r"^\|\s*\d+>\s?", "", l).strip() for l in t.splitlines()
                  if re.match(r"^\|\s*\d+>", l)]
            ec = [l for l in physics("\n".join(ec)).splitlines() if "END OF INPUT" not in l]
            fr, s = re.match(r"field_\d+_(\d+)_(\w+)\.out", p.name).groups()
            base = physics(read(INV / "sp_{}_{}.inp".format(fr, s))).splitlines()
            extra = [l for l in ec if l not in base]
            if not (all(l in ec for l in base) and len(extra) == 1 and extra[0].lower().startswith("%pointcharges")):
                echo_ok = False
        say("PASS" if echo_ok else "FAIL",
            "every s17 field input = its bare single-point input plus one %pointcharges line "
            "(run-control text set aside)")

    # ============================================================== 9 Phase 1 barriers
    head("9. PHASE 1 BARRIER PROVENANCE: ensemble_barriers.tsv vs path_summary.txt")
    worst, n9, bad9, below = 0.0, 0, [], []
    for r in tsv:
        f = r["frame"]
        p = ENS / "frame_{}".format(f) / "path_summary.txt"
        if not p.exists() or not r.get("barrier"):
            continue
        t = read(p)
        src = re.search(r"^# from .*/frame_(\d+)/neb\.out", t, re.M)
        if not src or src.group(1) != f:
            bad9.append(f + "(source path)")
        imgs = []
        for l in t.splitlines():
            m = re.match(r"^\s*(\d+)\s+([-\d.]+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", l)
            if m:
                imgs.append((int(m.group(1)), float(m.group(4)), "<= CI" in l))
        ci = [i for i in imgs if i[2]]
        if len(ci) != 1:
            bad9.append(f + "(CI rows {})".format(len(ci)))
            continue
        d = abs(ci[0][1] - float(r["barrier"]))
        worst, n9 = max(worst, d), n9 + 1
        if d > 0.005:
            bad9.append(f)
        pre = [i[1] for i in imgs if i[0] < ci[0][0]]
        if pre and min(pre) < -0.5:
            below.append("{}({:+.2f})".format(f, min(pre)))
    say("PASS" if not bad9 else "FAIL",
        "barrier = the '<= CI' dE of its own path_summary.txt: {} frames, max |diff| {:.3f} "
        "kcal/mol {}".format(n9, worst, bad9[:6] or ""))
    say("INFO", "{} of {} paths have an image before the CI more than 0.5 kcal/mol BELOW image 0, "
        "so the barrier is measured from a reactant end that is not the band's lowest point: {}".format(
            len(below), n9, " ".join(below)))

    # ============================================================== summary
    head("SUMMARY")
    print("  " + "   ".join("{} {}".format(k, counts[k]) for k in ("PASS", "FAIL", "GAP", "INFO")))
    return 1 if counts["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main())
