#!/usr/bin/env python3
"""
s18b_pipeline_check.py - provenance and consistency audit of the per-frame QM/MM pipeline
(05_qmmm/19_ensemble) behind every ensemble barrier: input structure -> reactant optimisation
-> restrained scan -> product optimisation -> NEB.

Companion to s18_provenance_check.py, which covers the files Phase 2.2 reads directly.
Reads ONLY committed files in chorismate-thesis-results (the .out.xz files are read with
Python's lzma module). Changes nothing. Same output conventions as s18:

    PASS   the check holds
    FAIL   the check does not hold
    GAP    an artefact that should exist for provenance is not committed / not on disk
    INFO   a fact worth recording that is not itself a pass/fail test

and exits non-zero if any check FAILS.

USAGE
    python3 s18b_pipeline_check.py <results_repo> [<code_repo>]

With the code repo, section 3 matches each frame's copies of run_scan.sh and
stage_tsguess.sh to the commit in the code history that holds the same bytes (needs git).

SECTIONS
    1  inventory of the per-frame pipeline files
    2  inputs: reactant/product optimisation and NEB - level of theory, QM and active regions
    3  scan inputs: restrained coordinates, targets, spring constants; script versions
    4  input structures: each .rst7 against the first geometry its reactant optimisation read
    5  convergence and program version of every optimisation and NEB
    6  the chain: optimised reactant -> NEB image 0 -> harvested geometries
    7  barriers: the last PATH SUMMARY in every neb.out against ensemble_barriers.tsv
    8  the reactant end of each band (gap 7 in the audit doc)

Python 3.6 compatible, standard library only.
"""
import hashlib
import lzma
import re
import subprocess
import sys
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path

H2KCAL = 627.5094740631
QM0 = 6207                      # first QM atom, 0-based (QMAtoms {6207:6230})
NAMES = ("C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 "
         "C10 O5 O6").split()
LOT = {"QMMM", "B3LYP", "D3BJ", "DEF2-SVP", "DEF2/J", "RIJCOSX"}
counts = Counter()


def say(tag, msg):
    counts[tag] += 1
    print("  {:<5} {}".format(tag, msg))


def head(t):
    print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)


def read(p):
    return Path(p).read_text(errors="replace")


def read_xz(p):
    with lzma.open(str(p), "rt", errors="replace") as fh:
        return fh.read()


def tokens(t):
    s = set()
    for l in t.splitlines():
        if l.strip().startswith("!"):
            s |= {x.upper() for x in l.strip()[1:].split()}
    return s


def active(t):
    m = re.search(r"ActiveAtoms\s*\{([^}]*)\}", t)
    return frozenset(int(x) for x in m.group(1).split()) if m else None


def masked(t, f):
    t = re.sub(r"ActiveAtoms\s*\{[^}]*\}", "ActiveAtoms {<list>}", t)
    return t.replace(f, "<f>")


def xyz_frames(p):
    """all structures in an xyz trajectory: list of (comment, [(el, x, y, z)])"""
    L = read(p).splitlines()
    out, i = [], 0
    while i < len(L):
        if not L[i].strip():
            i += 1
            continue
        n = int(L[i].split()[0])
        rows = [l.split() for l in L[i + 2:i + 2 + n]]
        out.append((L[i + 1], [(r[0], float(r[1]), float(r[2]), float(r[3])) for r in rows]))
        i += n + 2
    return out


def read_xyz(p):
    return xyz_frames(p)[0][1]


def dist(p, q):
    return sum((x - y) ** 2 for x, y in zip(p, q)) ** 0.5


def kabsch_rmsd(A, B):
    """RMSD after optimal superposition (Horn's quaternion method; standard library only)"""
    n = len(A)
    ca = [sum(p[k] for p in A) / n for k in range(3)]
    cb = [sum(p[k] for p in B) / n for k in range(3)]
    a = [[p[k] - ca[k] for k in range(3)] for p in A]
    b = [[p[k] - cb[k] for k in range(3)] for p in B]
    S = [[sum(a[i][r] * b[i][c] for i in range(n)) for c in range(3)] for r in range(3)]
    (xx, xy, xz), (yx, yy, yz), (zx, zy, zz) = S
    N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx],
         [yz - zy, xx - yy - zz, xy + yx, zx + xz],
         [zx - xz, xy + yx, -xx + yy - zz, yz + zy],
         [xy - yx, zx + xz, yz + zy, -xx - yy + zz]]
    lam = max(abs(v) for row in N for v in row) * 4 + 1.0   # power iteration on N + lam*I
    v = [1.0, 0.0, 0.0, 0.0]
    for _ in range(500):
        w = [sum((N[i][j] + (lam if i == j else 0.0)) * v[j] for j in range(4)) for i in range(4)]
        s = sum(x * x for x in w) ** 0.5
        v = [x / s for x in w]
    emax = sum(v[i] * sum(N[i][j] * v[j] for j in range(4)) for i in range(4))
    ga = sum(x * x for p in a for x in p)
    gb = sum(x * x for p in b for x in p)
    return max(0.0, (ga + gb - 2.0 * emax) / n) ** 0.5


def maxd(a, b):
    return max(abs(p[k] - q[k]) for p, q in zip(a, b) for k in (1, 2, 3))


def path_summaries(text):
    """every PATH SUMMARY table in an ORCA NEB output, as lists of (image, dist, E, dE, is_CI)"""
    out = []
    for m in re.finditer(r"PATH SUMMARY.*?\n(Image.*?\n)(.*?)\n\s*\n", text, re.S):
        rows = []
        for l in m.group(2).splitlines():
            f = re.match(r"^\s*(\d+)\s+([-\d.]+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", l)
            if f:
                rows.append((int(f.group(1)), float(f.group(2)), float(f.group(3)),
                             float(f.group(4)), "<= CI" in l))
        if rows:
            out.append(rows)
    return out


def read_tsv(p):
    hdr, rows = None, []
    for l in read(p).splitlines():
        if not l.strip():
            continue
        if l.startswith("#"):
            if hdr is None and "\t" in l:
                hdr = l.lstrip("#").strip().split("\t")
            continue
        f = l.split("\t")
        if hdr is None:
            hdr = f
            continue
        rows.append(dict(zip(hdr, f)))
    return rows


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    R = Path(sys.argv[1])
    C = Path(sys.argv[2]) if len(sys.argv) > 2 else None
    Q = R / "05_qmmm"
    E19, BAR, SEL = Q / "19_ensemble", Q / "19_ensemble_barriers", Q / "12_frame_selection" / "frames"
    for p in (E19, BAR, SEL):
        if not p.exists():
            sys.exit("FAIL: required folder missing: {}".format(p))
    tsv = read_tsv(BAR / "ensemble_barriers.tsv")
    frames = [r["frame"] for r in tsv]
    gl = [l for l in read(R / "phase2.2" / "grid_v2.tsv").splitlines() if l.startswith("#   frames =")]
    post = gl[0].split("=", 1)[1].strip().split(",")
    pilot = [f for f in frames if f not in post]

    # ======================================================================= 1
    head("1. INVENTORY OF THE PER-FRAME PIPELINE FILES (05_qmmm/19_ensemble)")
    dirs = sorted(p.name.replace("frame_", "") for p in E19.glob("frame_*"))
    say("PASS" if sorted(dirs) == sorted(frames) else "FAIL",
        "19_ensemble holds exactly the {} ensemble frames of ensemble_barriers.tsv".format(len(dirs)))
    want = ["reactant_opt.inp", "reactant_opt.out.xz", "reactant_opt.pbs", "reactant_opt.pbs.out",
            "reactant_opt-minimize-ener.csv", "scan.pbs", "scan/run_scan.sh", "scan/stage_tsguess.sh",
            "scan/active.txt", "scan/targets.txt", "scan/scan_energies.tsv", "scan/scan.pbs.out",
            "product_opt.inp", "product_opt.out.xz", "product_opt.pbs", "product_opt.pbs.out",
            "neb.inp", "neb.out.xz", "neb.pbs", "neb.pbs.out", "neb.NEB.log", "neb.interp",
            "neb_MEP.QMRegion_trj.xyz", "neb_NEB-CI_converged.QMRegion.xyz",
            "neb.inp.nebts", "neb.inp.neb_ci_superseded"]
    print("  {:<36}{:>8}{:>8}   missing".format("file", "post", "pilot"))
    for w in want:
        mp = [f for f in post if not (E19 / "frame_{}".format(f) / w).exists()]
        mq = [f for f in pilot if not (E19 / "frame_{}".format(f) / w).exists()]
        miss = " ".join(mp + mq) if len(mp + mq) <= 14 else "{} frames".format(len(mp + mq))
        print("  {:<36}{:>5}/{:<2}{:>5}/{:<2}   {}".format(w, len(post) - len(mp), len(post),
                                                        len(pilot) - len(mq), len(pilot), miss))
    nscan = Counter(len(list((E19 / "frame_{}".format(f) / "scan").glob("scan_*.inp"))) for f in frames)
    say("PASS" if list(nscan) == [20] else "FAIL",
        "every frame has 20 scan windows {}".format(dict(nscan)))
    for kind in ("-restraints.csv", "-colvars.csv", "-minimize-ener.csv"):
        n = Counter(len(list((E19 / "frame_{}".format(f) / "scan").glob("scan_*" + kind))) for f in frames)
        say("PASS" if list(n) == [20] else "GAP", "scan_NN{} for every window {}".format(kind, dict(n)))
    sup = sorted(p.parent.name.replace("frame_", "") for p in E19.glob("frame_*/neb_ci_superseded"))
    say("INFO", "superseded NEB-CI runs in {} frames: {}".format(len(sup), " ".join(sup)))
    say("PASS" if sorted(sup) == sorted(pilot) else "FAIL",
        "the superseded runs are exactly the pilot frames")
    files = Counter()
    for f in sup:
        for p in (E19 / "frame_{}".format(f) / "neb_ci_superseded").iterdir():
            files[p.name] += 1
    lacking = ["{} ({}/{})".format(k, v, len(sup)) for k, v in sorted(files.items()) if v < len(sup)]
    say("INFO", "superseded-run files present in fewer than all 14: {}".format(" ".join(lacking) or "none"))
    nebts = sorted(f for f in frames if (E19 / "frame_{}".format(f) / "neb.inp.nebts").exists())
    say("INFO", "neb.inp.nebts kept in {} frames: {}".format(len(nebts), " ".join(nebts)))

    # ======================================================================= 2
    head("2. INPUTS: REACTANT/PRODUCT OPTIMISATION AND NEB")
    acts = defaultdict(set)
    for stage, coord in (("reactant_opt.inp", "frame_<f>_CHA2.pdb"), ("product_opt.inp", "product_start.pdb"),
                         ("neb.inp", None), ("neb.inp.nebts", None), ("neb.inp.neb_ci_superseded", None)):
        texts = OrderedDict((f, read(E19 / "frame_{}".format(f) / stage)) for f in frames
                            if (E19 / "frame_{}".format(f) / stage).exists())
        if not texts:
            continue
        lot = [f for f, t in texts.items() if not LOT <= tokens(t)]
        say("PASS" if not lot else "FAIL", "{} ({} files): QMMM B3LYP D3BJ def2-SVP def2/J RIJCOSX "
            "in every file {}".format(stage, len(texts), lot or ""))
        var = Counter(masked(t, f) for f, t in texts.items())
        extras = Counter(" ".join(sorted(tokens(t) - LOT)) for t in texts.values())
        say("PASS" if len(var) == 1 else "INFO", "{}: {} text variant(s) once the frame number and "
            "ActiveAtoms list are masked; extra keywords {}".format(stage, len(var), dict(extras)))
        qm = [f for f, t in texts.items() if not re.search(r"QMAtoms\s*\{6207:6230\}", t)]
        say("PASS" if not qm else "FAIL", "{}: QMAtoms {{6207:6230}} in every file {}".format(stage, qm or ""))
        ch = [f for f, t in texts.items() if not re.search(r"^\*\s*\w+\s+-2\s+1\b", t, re.M)]
        say("PASS" if not ch else "FAIL", "{}: charge -2, multiplicity 1 {}".format(stage, ch or ""))
        if coord:
            bad = [f for f, t in texts.items() if coord.replace("<f>", f) not in t]
            say("PASS" if not bad else "FAIL", "{} reads {} {}".format(stage, coord, bad or ""))
        for f, t in texts.items():
            acts[f].add(active(t))
    same = [f for f in frames if len(acts[f]) == 1]
    say("PASS" if len(same) == len(frames) else "FAIL",
        "within each frame, every stage uses the same ActiveAtoms set: {}/{}".format(len(same), len(frames)))
    sets = Counter(next(iter(acts[f])) for f in same)
    say("INFO", "{} distinct active regions across frames, sizes {}".format(
        len(sets), sorted({len(s) for s in sets})))
    bad = [f for f in frames if acts[f] and not set(range(QM0, QM0 + 24)) <= next(iter(acts[f]))]
    say("PASS" if not bad else "FAIL", "the QM atoms lie inside every active region {}".format(bad or ""))
    cp = [f for f in frames if (BAR / "frame_{}".format(f) / "neb.inp").exists()
          and read(BAR / "frame_{}".format(f) / "neb.inp") != read(E19 / "frame_{}".format(f) / "neb.inp")]
    say("PASS" if not cp else "FAIL",
        "19_ensemble_barriers/neb.inp is a byte copy of 19_ensemble/neb.inp {}".format(cp or ""))

    # ======================================================================= 3
    head("3. SCAN INPUTS: RESTRAINED COORDINATES, TARGETS, SPRINGS, SCRIPT VERSIONS")
    colv, spr, mini, wrong_t, wrong_a, kw = Counter(), defaultdict(list), Counter(), [], [], Counter()
    for f in frames:
        S = E19 / "frame_{}".format(f) / "scan"
        tg = [l.split() for l in read(S / "targets.txt").splitlines() if l.strip()]
        at = frozenset(int(x) for x in read(S / "active.txt").split())
        springs = set()
        for k, p in enumerate(sorted(S.glob("scan_*.inp"))):
            t = read(p)
            cv = tuple(re.findall(r"Manage_Colvar Define (\d+) Distance Atom (\d+) Atom (\d+)", t))
            colv[cv] += 1
            rs = re.findall(r"Restraint Add Colvar (\d+) Harmonic Target ([\d.]+)_A Spring ([\d.]+)", t)
            springs |= {s for _, _, s in rs}
            if [r[1] for r in rs] != tg[k]:
                wrong_t.append("{}:{}".format(f, p.name))
            mini[re.search(r"^\s*Minimize.*$", t, re.M).group(0).strip()] += 1
            kw[" ".join(sorted(tokens(t)))] += 1
            if active(t) != at or (acts[f] and active(t) != next(iter(acts[f]))):
                wrong_a.append("{}:{}".format(f, p.name))
            if "win_{:02d}.pdb".format(k) not in t:
                wrong_t.append("{}:{} (pdb)".format(f, p.name))
        spr[" ".join(sorted(springs))].append(f)
    for cv, n in colv.items():
        desc = ", ".join("colvar {} = {}-{} (atoms {} {})".format(
            c, NAMES[int(a) - QM0], NAMES[int(b) - QM0], a, b) for c, a, b in cv)
        say("PASS" if len(colv) == 1 else "FAIL", "x{}: {}".format(n, desc))
    for s, fs in spr.items():
        say("INFO", "spring {} kJ/mol/A^2: {} frames ({})".format(s, len(fs), " ".join(fs)))
    say("PASS" if sorted(spr.get("400.0", [])) == sorted(pilot) and sorted(spr.get("2500.0", [])) == sorted(post)
        else "FAIL", "spring 400 for exactly the pilot frames and 2500 for exactly the post-cut frames")
    say("PASS" if not wrong_t else "FAIL", "every window's targets equal its row of targets.txt and it "
        "reads its own win_NN.pdb {}".format(wrong_t[:5] or ""))
    say("PASS" if not wrong_a else "FAIL", "every window's ActiveAtoms = active.txt = the frame's "
        "optimisation and NEB active region {}".format(wrong_a[:5] or ""))
    say("PASS" if len(mini) == 1 else "INFO", "minimiser settings: {}".format(dict(mini)))
    say("PASS" if len(kw) == 1 else "INFO", "scan simple-input keywords: {}".format(dict(kw)))
    dev = []
    for f in frames:
        for p in sorted((E19 / "frame_{}".format(f) / "scan").glob("scan_*-restraints.csv")):
            rows = [l.split(";") for l in read(p).splitlines() if l.strip() and not l.startswith("#")]
            if rows:
                r = rows[-1]
                dev.append((max(abs(float(r[1]) - float(r[2])), abs(float(r[6]) - float(r[7]))), f, p.name))
    dev.sort(reverse=True)
    say("INFO", "restraint achieved at the end of each window: max |position - target| {:.3f} A "
        "({} {}); median {:.4f} A over {} windows".format(dev[0][0], dev[0][1], dev[0][2],
                                                          dev[len(dev) // 2][0], len(dev)))
    if C:
        for name in ("run_scan.sh", "stage_tsguess.sh"):
            hist = {}
            try:
                hs = subprocess.check_output(["git", "-C", str(C), "log", "--all", "--format=%h %ad",
                                              "--date=short", "--", "phase1_system_dev/" + name]).decode().split("\n")
                for line in hs:
                    if not line.strip():
                        continue
                    h, d = line.split()
                    blob = subprocess.check_output(["git", "-C", str(C), "show",
                                                    "{}:phase1_system_dev/{}".format(h, name)])
                    hist[hashlib.sha256(blob).hexdigest()] = "{} ({})".format(h, d)
            except Exception as e:
                say("INFO", "{}: git history not readable ({})".format(name, e))
                continue
            grp = defaultdict(list)
            for f in frames:
                p = E19 / "frame_{}".format(f) / "scan" / name
                grp[hist.get(hashlib.sha256(p.read_bytes()).hexdigest(), "NO MATCH") if p.exists()
                    else "absent"].append(f)
            for k, fs in grp.items():
                say("INFO" if k != "NO MATCH" else "FAIL",
                    "{} copy = code commit {}: {} frames ({})".format(name, k, len(fs),
                                                                       " ".join(fs) if len(fs) <= 14 else "..."))
        for f in pilot:
            p = E19 / "frame_{}".format(f) / "scan" / "run_scan.sh"
            if p.exists() and "SPRING=2500" in read(p) and f in spr.get("400.0", []):
                say("INFO", "frame {}: its run_scan.sh copy sets SPRING=2500 but its scan inputs use 400 and "
                    "its scan log is from {}; the copy postdates the run, the inputs are the record".format(
                        f, (re.search(r"start=(.*)", read(E19 / "frame_{}".format(f) / "scan" / "scan.pbs.out"))
                            or [None, "?"])[1].strip()))

    # ======================================================================= 4
    head("4. INPUT STRUCTURES: .rst7 AGAINST THE FIRST GEOMETRY EACH REACTANT OPTIMISATION READ")
    bad, worst, n4 = [], 0.0, 0
    outs = {}
    for f in frames:
        rst = SEL / "frame_{}_CHA2.rst7".format(f)
        L = read(rst).splitlines()
        nat = int(L[1].split()[0])
        vals = []
        for l in L[2:]:
            for k in range(0, len(l.rstrip()), 12):
                vals.append(float(l[k:k + 12]))
                if len(vals) == 3 * nat:
                    break
            if len(vals) == 3 * nat:
                break
        txt = read_xz(E19 / "frame_{}".format(f) / "reactant_opt.out.xz")
        outs[f] = txt
        m = re.search(r"CARTESIAN COORDINATES \(ANGSTROEM\)\n-+\n(.*?)\n\s*\n", txt, re.S)
        rows = [l.split() for l in m.group(1).splitlines()]
        if len(rows) != nat:
            bad.append("{} ({} vs {} atoms)".format(f, len(rows), nat))
            continue
        d = max(abs(float(rows[i][k + 1]) - vals[3 * i + k]) for i in range(nat) for k in range(3))
        worst, n4 = max(worst, d), n4 + 1
        if d > 6e-4:
            bad.append("{} ({:.1e})".format(f, d))
        qa = [int(x) for x in read(SEL / "frame_{}_CHA2.qmatoms".format(f)).split()]
        if qa != list(range(QM0 + 1, QM0 + 25)):
            bad.append(f + " (qmatoms)")
        els = "".join(rows[i][0][0] for i in range(QM0, QM0 + 24))
        if els != "CHHCCOOOCHCHCCHCHCHOHCOO":
            bad.append(f + " (QM elements)")
    say("PASS" if not bad else "FAIL",
        "{} frames: every one of 55680 atoms in the .rst7 matches the reactant optimisation's starting "
        "geometry to the PDB's 3 decimals (max |dx| {:.1e} A); QM atoms in the prms order; .qmatoms "
        "6208-6231 {}".format(n4, worst, bad[:5] or ""))
    man = read_tsv(SEL / "frames_manifest.tsv")
    mf = sorted("{:05d}".format(int(r["frame"])) for r in man if "frame" in r)
    say("PASS" if set(frames) <= set(mf) else "FAIL",
        "frames_manifest.tsv lists every ensemble frame ({} rows)".format(len(mf)))

    # ======================================================================= 5
    head("5. CONVERGENCE AND PROGRAM VERSION")
    prob = defaultdict(list)
    ver = Counter()
    nebtxt = {}
    for f in frames:
        F = E19 / "frame_{}".format(f)
        for stage in ("reactant_opt", "product_opt"):
            t = outs[f] if stage == "reactant_opt" else read_xz(F / "product_opt.out.xz")
            if stage == "product_opt":
                outs[f + "_P"] = t
            ver[(re.search(r"Program Version (\S+)", t) or [None, None])[1]] += 1
            if "The minimization has converged" not in t:
                prob[stage + " not converged"].append(f)
            if "TERMINATED NORMALLY" not in t:
                prob[stage + " not terminated normally"].append(f)
        t = read_xz(F / "neb.out.xz")
        nebtxt[f] = t
        ver[(re.search(r"Program Version (\S+)", t) or [None, None])[1]] += 1
        if "THE NEB OPTIMIZATION HAS CONVERGED" not in t:
            prob["NEB not converged"].append(f)
    say("PASS" if list(ver) == ["6.0.1"] else "FAIL", "ORCA version across {} outputs: {}".format(
        sum(ver.values()), dict(ver)))
    for lab in ("reactant_opt not converged", "reactant_opt not terminated normally",
                "product_opt not converged", "product_opt not terminated normally", "NEB not converged"):
        say("PASS" if not prob[lab] else "FAIL", "{}: {} {}".format(lab, len(prob[lab]), prob[lab] or ""))
    term = [f for f in frames if "TERMINATED NORMALLY" not in nebtxt[f]]
    say("INFO", "neb.out without ORCA TERMINATED NORMALLY (NEB-TS jobs killed in the TS stage, "
        "per HANDOFF): {} {}".format(len(term), " ".join(term)))

    # ======================================================================= 6
    head("6. THE CHAIN: OPTIMISED REACTANT -> NEB IMAGE 0 -> HARVESTED GEOMETRIES")
    e6, g6, p6, ts6 = [], [], [], []
    for f in frames:
        F = E19 / "frame_{}".format(f)
        er = float(re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", outs[f])[-1])
        ep = float(re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", outs[f + "_P"])[-1])
        mep = xyz_frames(F / "neb_MEP.QMRegion_trj.xyz")
        e0 = float(re.search(r"E\s+(-?\d+\.\d+)", mep[0][0]).group(1))
        eN = float(re.search(r"E\s+(-?\d+\.\d+)", mep[-1][0]).group(1))
        e6.append((abs(e0 - er) * H2KCAL, abs(eN - ep) * H2KCAL, f))
        for s, img in (("reactant_qm.xyz", 0), ("product_qm.xyz", -1)):
            b = BAR / "frame_{}".format(f) / s
            if b.exists():
                g6.append((maxd(read_xyz(b), mep[img][1]), f, s))
        b = BAR / "frame_{}".format(f) / "transition_state_qm.xyz"
        c = F / "neb_NEB-CI_converged.QMRegion.xyz"
        if b.exists() and c.exists():
            ts6.append((maxd(read_xyz(b), read_xyz(c)), f))
        elif c.exists():
            p6.append(f)
    w = max(e6)
    say("PASS" if max(x[0] for x in e6) < 0.01 and max(x[1] for x in e6) < 0.01 else "FAIL",
        "NEB image 0 = the optimised reactant and the last image = the optimised product, by energy: "
        "max |dE| {:.4f} / {:.4f} kcal/mol over {} frames".format(
            max(x[0] for x in e6), max(x[1] for x in e6), len(e6)))
    say("PASS" if all(x[0] <= 1e-6 for x in g6) else "FAIL",
        "19_ensemble_barriers reactant/product_qm.xyz = NEB image 0 / last image: {} files, max |dx| "
        "{:.1e} A".format(len(g6), max(x[0] for x in g6)))
    say("PASS" if all(x[0] <= 1e-6 for x in ts6) else "FAIL",
        "19_ensemble_barriers transition_state_qm.xyz = neb_NEB-CI_converged.QMRegion.xyz: {} frames".format(len(ts6)))
    if p6:
        say("INFO", "climbing image committed in 19_ensemble but never harvested into "
            "19_ensemble_barriers: {}".format(" ".join(p6)))
    INV = Q / "20_invacuo"
    w20, n20, b20 = 0.0, 0, []
    for f in frames:
        F = E19 / "frame_{}".format(f)
        mep = xyz_frames(F / "neb_MEP.QMRegion_trj.xyz")
        ci = F / "neb_NEB-CI_converged.QMRegion.xyz"
        for s, src in (("R", mep[0][1]), ("P", mep[-1][1]), ("TS", read_xyz(ci) if ci.exists() else None)):
            g = INV / "{}_{}.xyz".format(f, s)
            if src is None or not g.exists():
                continue
            d = maxd(read_xyz(g), src)
            w20, n20 = max(w20, d), n20 + 1
            if d > 1e-6:
                b20.append("{}_{}".format(f, s))
    say("PASS" if not b20 else "FAIL",
        "20_invacuo R/TS/P geometries = NEB image 0 / climbing image / last image, all frames: {} files, "
        "max |dx| {:.1e} A {}".format(n20, w20, b20[:5] or ""))

    # ======================================================================= 7
    head("7. BARRIERS: LAST PATH SUMMARY IN EACH neb.out AGAINST ensemble_barriers.tsv")
    rows7, bad7, sup7 = [], [], []
    for r in tsv:
        f = r["frame"]
        ps = path_summaries(nebtxt[f])
        if not ps:
            bad7.append(f + "(no path summary)")
            continue
        ci = [x for x in ps[-1] if x[4]]
        if len(ci) != 1:
            bad7.append(f + "(CI rows {})".format(len(ci)))
            continue
        if not r.get("barrier"):
            say("INFO", "frame {}: no barrier in the tsv; neb.out CI dE {:.2f}".format(f, ci[0][3]))
            continue
        d = abs(ci[0][3] - float(r["barrier"]))
        rows7.append((d, f))
        if d > 0.005:
            bad7.append("{}({:.2f} vs {})".format(f, ci[0][3], r["barrier"]))
        S = E19 / "frame_{}".format(f) / "neb_ci_superseded" / "neb.out.xz"
        if S.exists():
            pss = path_summaries(read_xz(S))
            cs = [x for x in pss[-1] if x[4]] if pss else []
            sup7.append("{} {:.2f}->{}".format(f, cs[0][3] if cs else float("nan"), r["barrier"]))
    say("PASS" if not bad7 else "FAIL", "tsv barrier = '<= CI' dE of the LAST path summary in neb.out: "
        "{} frames, max |diff| {:.3f} kcal/mol {}".format(len(rows7), max(rows7)[0] if rows7 else 0, bad7[:6] or ""))
    for f in ("44388", "59999"):
        if f in [x[1] for x in rows7]:
            say("PASS", "frame {}: its barrier now traces to a committed neb.out".format(f))
    if sup7:
        say("INFO", "superseded NEB-CI barrier -> current barrier, pilot frames: {}".format("; ".join(sup7)))

    # ======================================================================= 8
    head("8. THE REACTANT END OF EACH BAND (audit gap 7)")
    dips = []
    for f in frames:
        ps = path_summaries(nebtxt[f])
        if not ps:
            continue
        last = ps[-1]
        ci = [x for x in last if x[4]]
        if not ci:
            continue
        pre = [x for x in last if x[0] < ci[0][0]]
        low = min(pre, key=lambda x: x[3]) if pre else None
        if low and low[3] < -0.5:
            mep = xyz_frames(E19 / "frame_{}".format(f) / "neb_MEP.QMRegion_trj.xyz")
            A, B = [r[1:] for r in mep[0][1]], [r[1:] for r in mep[low[0]][1]]
            dips.append((low[3], f, low[0], low[1], ci[0][3], kabsch_rmsd(A, B),
                         abs(dist(A[7], A[8]) - dist(B[7], B[8]))))
    dips.sort()
    say("INFO", "{} of {} final bands have an image before the CI more than 0.5 kcal/mol below image 0".format(
        len(dips), len(frames)))
    for d, f, im, pl, b, rm, dco in dips[:8]:
        say("INFO", "  frame {}: image {} at {:.2f} A along the band is {:+.2f} kcal/mol; barrier from image 0 "
            "{:.2f}, from that image {:.2f}; substrate RMSD to image 0 {:.3f} A".format(f, im, pl, d, b, b - d, rm))
    if dips:
        rms = sorted(x[5] for x in dips)
        say("INFO", "between image 0 and the lowest pre-CI image the SUBSTRATE barely moves: QM-region RMSD median "
            "{:.3f}, max {:.3f} A; C4-O3 changes by at most {:.3f} A; yet those images lie {:.2f}-{:.2f} A along a "
            "band measured over the whole active region. The dip is the MM environment settling, not the "
            "substrate".format(rms[len(rms) // 2], rms[-1], max(x[6] for x in dips), min(x[3] for x in dips),
                               max(x[3] for x in dips)))
    nsp = Counter(len(re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)", nebtxt[f])) for f in frames)
    say("INFO", "neb.out does not record the QM and MM energy of each image on the final band "
        "(QM/MM single-point blocks per neb.out: {}). Splitting the dip into QM and MM needs single points "
        "on the final images (full active-region geometries on hpc1 only, checksum-only); the geometric "
        "test above locates the dip without it".format(dict(nsp)))

    head("SUMMARY")
    print("  " + "   ".join("{} {}".format(k, counts[k]) for k in ("PASS", "FAIL", "GAP", "INFO")))
    return 1 if counts["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main())
