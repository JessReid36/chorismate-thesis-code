#!/usr/bin/env python3
"""
cs5_align.py - Stage 5 (phase2.2, STAGE5_CRITERIA.txt): the production geometries for the grid, L and the A matrix.
The substrate (24 atoms) at the SP package's TS, R_IRC and P of every complete frame, taken from the package's own
final geometries, then superimposed as s16_align_frames.py does.

  python3 cs5_align.py extract SD TABLE REPORT OUT_EXTRACTED   (on hpc1, or on a copy of the run folders; standard library)
  python3 cs5_align.py align   OUT_EXTRACTED OUT_ALIGNED  (numpy)

EXTRACT
  SD     the folder holding the package's run folders (f<FRAME>_r8_a1; SP_AMENDMENT2.txt)
  TABLE  spa_package.py's per-frame table, written by its FINAL run (first header line ends "FINAL"; A2.9)
  REPORT that run's printed report, saved: it must say "RESULT: FINAL (nothing undecided)" and "THE PHASE 1 BARRIER
         (A2.3; audit item C3): the conditions are met" (STAGE5_CRITERIA.txt, WHEN IT RUNS); its sha256 is recorded
  Frames: every row with surroundings "md" and state "complete" (SP_PLAN.txt QUANTITIES: TS, TS frequencies, IRC,
  R_IRC and R_IRC frequencies passed; SP_AMENDMENT2.txt A2.2). The surroundings-check runs (S = ci) are not used.
  Per frame, from the run folder:
    TS  ts/final.xyz      the converged alternation's last core geometry (spa_stage.py), whose sha256 must equal
                          ts/STAGE_RESULT.txt's xyz_sha256 with status=PASS
    R   r_irc/final.xyz   likewise (r_irc/STAGE_RESULT.txt), the reactant on the TS's own path (B_sp's reactant)
    P   p_irc/final.xyz   if p_irc/STAGE_RESULT.txt has status=PASS (sha256 checked); otherwise the IRC's product end
                          irc/p_irc_start.xyz, whose sha256 must equal irc/STAGE_RESULT.txt's p_irc_start_sha256 and
                          whose end check_irc labelled "product". P enters only the grid's surface (the union of R, TS
                          and P) and the optional V_P column; neither the a-matrix nor L uses it.
  The substrate is atoms 6207-6230 (0-based, ORCA's QMAtoms {6207:6230}) of each full-system file, copied as text
  (full precision). Checked: 24 atoms with the substrate's element order; every full-system file has the same atom
  count; the breaking O3-C4 and forming C1-C6 distances equal the stage record's break_A and form_A (3 decimals), and
  for P from the IRC end, check_irc's recorded end (3 decimals).
  Writes OUT_EXTRACTED/frame_<F>/{R,TS,P}.xyz and EXTRACT.txt (per frame: folder, sources, sha256, P's source).

ALIGN (as s16_align_frames.py, whose kabsch() is reproduced below unchanged)
  One rotation and translation per frame, fitted on the rigid core (C2, O3 and the six ring carbons: indices
  3, 7, 8, 10, 12, 13, 15, 17) of the frame's R onto the reference frame's R, and applied unchanged to its TS and P, so
  the reactant-to-TS motion is kept. Reference: the lowest-numbered frame (frame 20000, as D3, when it is complete).
  Writes OUT_ALIGNED/frame_<F>/{reactant_qm,transition_state_qm,product_qm}.xyz (the names cs5_grid.py and s15
  read), transform.txt (aligned = original @ R + t) and ALIGNMENT.txt (RMSD on the core before and after, per frame).
  Checked: the rotation is proper and orthonormal, and the inverse transform returns every original coordinate within
  1e-8 A.
Both modes refuse to write into an existing folder.
"""
import hashlib, os, sys

NAMES = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
QM0 = 6207                                   # ORCA's QMAtoms {6207:6230}, 0-based (the package's core.inp)
BREAK = (NAMES.index("O3"), NAMES.index("C4"))   # 6214-6215 in spa_stage.py
FORM = (NAMES.index("C1"), NAMES.index("C6"))    # 6207-6219
CORE = [3, 7, 8, 10, 12, 13, 15, 17]         # s16_align_frames.py CORE
STATE_OUT = {"R": "reactant_qm.xyz", "TS": "transition_state_qm.xyz", "P": "product_qm.xyz"}

def stop(m): sys.stderr.write("STOP: " + m + "\n"); sys.exit(1)
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
def record(p):
    if not os.path.isfile(p): stop("%s is missing" % p)
    return dict(l.split("=", 1) for l in text(p).splitlines() if "=" in l)
def dist(a, b): return sum((a[k] - b[k]) ** 2 for k in range(3)) ** 0.5

def substrate(p):
    """the 24 substrate lines of a full-system xyz, as text, with their coordinates and the file's atom count"""
    with open(p, encoding="utf-8", errors="replace") as fh:
        n = int(fh.readline().split()[0]); fh.readline()
        lines = []
        for i in range(n):
            l = fh.readline()
            if not l: stop("%s ends before its %d atoms" % (p, n))
            if QM0 <= i < QM0 + 24: lines.append(l.rstrip("\n"))
    if len(lines) != 24: stop("%s: fewer than %d atoms" % (p, QM0 + 24))
    el = [l.split()[0] for l in lines]
    if [e[0] for e in el] != [x[0] for x in NAMES] or any(len(e) != 1 for e in el):
        stop("%s: atoms %d-%d are not the substrate (%s)" % (p, QM0, QM0 + 23, " ".join(el)))
    X = [[float(v) for v in l.split()[1:4]] for l in lines]
    return n, lines, X

def text(p): return open(p, encoding="utf-8", errors="replace").read()

def extract(sd, table, report, out):
    if os.path.exists(out): stop("%s exists" % out)
    rep = text(report).splitlines()
    if "RESULT: FINAL (nothing undecided)" not in rep or not any(
            l.startswith("THE PHASE 1 BARRIER (A2.3; audit item C3): the conditions are met") for l in rep):
        stop("%s does not state a FINAL result with the Phase 1 barrier conditions met: Stage 5 waits (STAGE5_CRITERIA.txt)" % report)
    T = text(table).splitlines()
    if not T or not T[0].startswith("# spa_package.py") or not T[0].rstrip().endswith("FINAL"):
        stop("%s is not spa_package.py's FINAL table" % table)
    h = T[1].split("\t")
    for c in ("frame", "surroundings", "folder", "state"):
        if c not in h: stop("%s: no column %s" % (table, c))
    rows = [dict(zip(h, l.split("\t"))) for l in T[2:] if l.strip()]
    sel = sorted((r for r in rows if r["surroundings"] == "md" and r["state"] == "complete"), key=lambda r: int(r["frame"]))
    if not sel: stop("no complete frame (surroundings md) in %s" % table)
    if len(set(r["frame"] for r in sel)) != len(sel): stop("a frame appears twice")
    os.makedirs(out)
    log = ["# cs5_align.py extract (Stage 5, STAGE5_CRITERIA.txt); cs5_align.py sha256 %s" % sha(__file__),
           "# SD %s; table %s sha256 %s; report %s sha256 %s" % (sd, table, sha(table), report, sha(report)),
           "# frames (%d): %s" % (len(sel), " ".join(r["frame"] for r in sel)),
           "# substrate = atoms %d-%d (0-based) of each full-system file, copied as text" % (QM0, QM0 + 23),
           "frame\tfolder\tstate\tsource\tsha256 (full system)\tstage record\tbreak_A\tform_A"]
    natoms = None
    for r in sel:
        rd = os.path.join(sd, r["folder"])
        if r["folder"] != "f%s_r8_a1" % r["frame"]: stop("frame %s: folder %s is not f%s_r8_a1" % (r["frame"], r["folder"], r["frame"]))
        src = {}
        for s, step in (("TS", "ts"), ("R", "r_irc")):
            rec = record(os.path.join(rd, step, "STAGE_RESULT.txt"))
            if rec.get("status") != "PASS": stop("%s: %s has not passed (status %s)" % (rd, step, rec.get("status")))
            src[s] = (os.path.join(rd, step, "final.xyz"), rec, rec.get("xyz_sha256"), "%s/STAGE_RESULT.txt" % step)
        pr = os.path.join(rd, "p_irc", "STAGE_RESULT.txt")
        rec = record(pr) if os.path.isfile(pr) else {}
        if rec.get("status") == "PASS":
            src["P"] = (os.path.join(rd, "p_irc", "final.xyz"), rec, rec.get("xyz_sha256"), "p_irc/STAGE_RESULT.txt")
        else:
            irc = record(os.path.join(rd, "irc", "STAGE_RESULT.txt"))
            if irc.get("status") != "PASS": stop("%s: the IRC has not passed, yet the frame is complete" % rd)
            end = [k for k in ("end_F", "end_B") if irc.get(k, "").startswith("product")]
            if len(end) != 1 or irc.get("p_irc_from") != "irc_IRC_%s.xyz" % end[0][-1]:
                stop("%s: irc/STAGE_RESULT.txt does not name one product end as p_irc's source" % rd)
            b, f = irc[end[0]].split()[2], irc[end[0]].split()[4]
            src["P"] = (os.path.join(rd, "irc", "p_irc_start.xyz"), {"break_A": b, "form_A": f}, irc.get("p_irc_start_sha256"),
                        "irc/STAGE_RESULT.txt %s (p_irc %s)" % (end[0], rec.get("status", "not run")))
        d = os.path.join(out, "frame_" + r["frame"]); os.makedirs(d)
        for s in ("R", "TS", "P"):
            p, rec, want, what = src[s]
            if not want or sha(p) != want: stop("%s: sha256 differs from %s's record" % (p, what))
            n, lines, X = substrate(p)
            if natoms is None: natoms = n
            if n != natoms: stop("%s has %d atoms, not %d like the first file" % (p, n, natoms))
            br, fo = dist(X[BREAK[0]], X[BREAK[1]]), dist(X[FORM[0]], X[FORM[1]])
            if abs(br - float(rec["break_A"])) > 0.0006 or abs(fo - float(rec["form_A"])) > 0.0006:
                stop("%s: O3-C4 %.4f and C1-C6 %.4f differ from the record's %s and %s" % (p, br, fo, rec["break_A"], rec["form_A"]))
            with open(os.path.join(d, s + ".xyz"), "w") as fh:
                fh.write("24\nframe %s %s: substrate (atoms %d-%d) of %s/%s\n" % (r["frame"], s, QM0, QM0 + 23, r["folder"],
                                                                                 os.path.relpath(p, rd)) + "\n".join(lines) + "\n")
            log.append("\t".join([r["frame"], r["folder"], s, os.path.relpath(p, rd), want, what, "%.3f" % br, "%.3f" % fo]))
    open(os.path.join(out, "EXTRACT.txt"), "w").write("\n".join(log) + "\n")
    print("extracted %d frames x R, TS, P into %s (full systems of %d atoms); P from the IRC end in %d frame(s)" % (
        len(sel), out, natoms, sum(1 for l in log if "\tP\tirc/" in l)))

def kabsch(P, Q, idx):
    """s16_align_frames.py's kabsch(), unchanged: (R, t) such that P @ R + t is the aligned structure."""
    import numpy as np
    Pc, Qc = P[idx].mean(0), Q[idx].mean(0)
    A = (P[idx] - Pc).T @ (Q[idx] - Qc)
    V, S, W = np.linalg.svd(A)
    d = np.sign(np.linalg.det(V @ W))
    R = V @ np.diag([1.0, 1.0, d]) @ W
    return R, Qc - Pc @ R

def read24(p):
    import numpy as np
    L = open(p).read().splitlines()
    if int(L[0]) != 24: stop("%s: not 24 atoms" % p)
    return [l.split()[0] for l in L[2:26]], np.array([[float(v) for v in l.split()[1:4]] for l in L[2:26]])

def align(src, out):
    import numpy as np
    if os.path.exists(out): stop("%s exists" % out)
    frames = sorted((d[6:] for d in os.listdir(src) if d.startswith("frame_")), key=int)
    if not frames: stop("no frame_* in %s" % src)
    ref = frames[0]
    _, Xref = read24(os.path.join(src, "frame_" + ref, "R.xyz"))
    os.makedirs(out)
    rows = []
    for f in frames:
        G = {s: read24(os.path.join(src, "frame_" + f, s + ".xyz")) for s in ("R", "TS", "P")}
        R, t = kabsch(G["R"][1], Xref, CORE)
        if abs(np.linalg.det(R) - 1) > 1e-9 or np.abs(R @ R.T - np.eye(3)).max() > 1e-9: stop("frame %s: rotation not proper" % f)
        d = os.path.join(out, "frame_" + f); os.makedirs(d)
        np.savetxt(os.path.join(d, "transform.txt"), np.vstack([R, t.reshape(1, 3)]), fmt="%.12f",
                   header="rows 0-2: rotation R; row 3: translation t\naligned = original @ R + t ; original = (aligned - t) @ R.T")
        M = np.loadtxt(os.path.join(d, "transform.txt")); Rr, tr = M[:3], M[3]       # as written, as s15/cs5_amatrix read it
        for s in ("R", "TS", "P"):
            el, X = G[s]
            Y = X @ Rr + tr
            with open(os.path.join(d, STATE_OUT[s]), "w") as fh:
                fh.write("24\nframe %s %s, aligned onto frame %s R on the rigid core (%d atoms)\n" % (f, s, ref, len(CORE)))
                for e, (x, y, z) in zip(el, Y): fh.write("%-3s%18.10f%18.10f%18.10f\n" % (e, x, y, z))
            _, Yr = read24(os.path.join(d, STATE_OUT[s]))
            if np.abs((Yr - tr) @ Rr.T - X).max() > 1e-8: stop("frame %s %s: the inverse transform does not return the original" % (f, s))
        pre = np.sqrt(((G["R"][1][CORE] - Xref[CORE]) ** 2).sum(1).mean())
        post = np.sqrt(((G["R"][1][CORE] @ Rr + tr - Xref[CORE]) ** 2).sum(1).mean())
        rows.append((f, pre, post))
    post = [r[2] for r in rows]
    txt = ["cs5_align.py align (Stage 5, STAGE5_CRITERIA.txt); cs5_align.py sha256 %s" % sha(__file__),
           "source %s (EXTRACT.txt sha256 %s)" % (src, sha(os.path.join(src, "EXTRACT.txt")) if os.path.isfile(os.path.join(src, "EXTRACT.txt")) else "-"),
           "reference: frame %s, R (R_IRC)" % ref,
           "alignment set: rigid core, atom indices %s (C2, O3 and the six ring carbons; s16_align_frames.py)" % CORE,
           "one transform per frame, fitted on R and applied unchanged to TS and P",
           "frames (%d): %s" % (len(frames), " ".join(frames)),
           "RMSD on the core after alignment: min %.3f mean %.3f max %.3f A" % (min(post), sum(post) / len(post), max(post)),
           "", "frame\tRMSD_before_A\tRMSD_after_A"] + ["%s\t%.4f\t%.4f" % r for r in rows]
    open(os.path.join(out, "ALIGNMENT.txt"), "w").write("\n".join(txt) + "\n")
    print("aligned %d frames onto frame %s R; core RMSD after %.3f-%.3f A" % (len(frames), ref, min(post), max(post)))

if __name__ == "__main__":
    a = sys.argv[1:]
    if a[:1] == ["extract"] and len(a) == 5: extract(*a[1:])
    elif a[:1] == ["align"] and len(a) == 3: align(*a[1:])
    else: stop(__doc__.split("\n\n")[1])
