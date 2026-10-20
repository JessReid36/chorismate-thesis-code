#!/usr/bin/env python3
"""
neb_active_regions.py - the optimised enzyme around the substrate at the reactant and the transition state of each
Phase 1 band frame, from the NEB's own active-region files.

WHY
The site-effect map draws the substrate of frame 20000 as Phase 1 harvested it: the reactant is NEB image 0
(05_qmmm/19_ensemble_barriers/frame_N/reactant_qm.xyz = neb_MEP.QMRegion_trj.xyz image 0, s18b_pipeline_check.py
section 6) and the transition state is the climbing image (transition_state_qm.xyz = neb_NEB-CI_converged.QMRegion.xyz).
ORCA wrote the optimised region around it (QM atoms + active MM atoms, 1959 atoms) for the same two structures:
neb_MEP.activeRegion_trj.xyz (image 0) and neb_NEB-CI_converged.activeRegion.xyz. These were checksum-recorded, not
committed (05_qmmm/19_ensemble/CHECKSUMS_not_committed.sha256). Taking residues from them, mapped with the frame's
committed transform (phase2.2/aligned/frame_N/transform.txt), extracts them exactly as the substrate was extracted.

MODE
  python3 neb_active_regions.py check <results_repo> <folder holding 19_ensemble/frame_N/...> <report>
    For each of the 30 band frames (phase2.2/cs_stage5/exploratory/dev_A/A_dev.tsv): the active atoms of the frame's
    neb.inp; 10 images of 1959 atoms in the MEP file, 1959 in the CI file; every element the prmtop's; the 24 QM atoms
    of image 0, the last image and the CI equal to the committed reactant_qm.xyz, product_qm.xyz and
    transition_state_qm.xyz (1e-5 A); at frame 41786, every active MM atom equal to Stage 4d's env_R (pc/ENZ.pc.xz,
    from reactant.pdb: 1e-3 A, the PDB's rounding) and env_TS (pc/ENZTS.pc.xz, from the CI: 1e-5 A). Writes the
    report and exits 1 at the first failure. Reports, per frame, how far the six charged residues Stage 4d analyses
    move between the reactant and the transition state.
Python 3.6+, numpy.
"""
import lzma, os, re, sys
import numpy as np

QM0, NQM = 6207, 24
RES = [("Arg90", 217, "CZ"), ("Arg7", 134, "CZ"), ("Glu78", 205, "CD"), ("Arg116", 243, "CZ"),
       ("Arg63'", 63, "CZ"), ("Lys60'", 60, "NZ")]          # Stage 4d's six (cs4d_prepare.py), charge centres (cs4d_explore.py)
SYM = {1: "H", 6: "C", 7: "N", 8: "O", 11: "NA", 16: "S", 17: "CL"}
out = []


def say(s=""):
    print(s); out.append(s)


def stop(rep, s):
    say("STOP: " + s)
    open(rep, "w").write("\n".join(out) + "\n")
    sys.exit(1)


def prmtop(path):
    T = [l.rstrip() for l in open(path)]

    def flag(name):
        i = T.index("%FLAG " + name); fmt = T[i + 1]; j = i + 2; v = []
        w = int(fmt.split("a")[1].split(")")[0]) if "a" in fmt else int(fmt.split("I")[1].split(")")[0])
        while not T[j].startswith("%FLAG"):
            v += [T[j][k:k + w].strip() for k in range(0, len(T[j]), w)]; j += 1
        return [x for x in v if x]
    an, z = flag("ATOM_NAME"), [int(x) for x in flag("ATOMIC_NUMBER")]
    rp = [int(x) - 1 for x in flag("RESIDUE_POINTER")] + [len(an)]
    return an, z, rp


def frames_xyz(path):
    """every structure in an xyz file (one or many): list of (comment, elements, coordinates)"""
    L = open(path).read().splitlines(); k = 0; F = []
    while k < len(L) and L[k].strip():
        n = int(L[k].split()[0]); c = L[k + 1]
        rows = [l.split() for l in L[k + 2:k + 2 + n]]
        if len(rows) != n or any(len(r) < 4 for r in rows):
            raise ValueError("%s: structure %d is cut short" % (path, len(F) + 1))
        F.append((c, [r[0] for r in rows], np.array([[float(x) for x in r[1:4]] for r in rows])))
        k += 2 + n
    return F


def check(res, dl, rep):
    fr = open(res + "/phase2.2/cs_stage5/exploratory/dev_A/A_dev.tsv").read().split("# frames (30): ")[1].split("\n")[0].split()
    if len(fr) != 30: stop(rep, "A_dev.tsv does not list 30 frames")
    an, z, rp = prmtop(res + "/03_amber/tleap_build/complex_solvated.prmtop")
    say("NEB active regions of the 30 Phase 1 band frames, checked by neb_active_regions.py")
    move = []
    for f in fr:
        E = res + "/05_qmmm/19_ensemble/frame_%s/" % f; B = res + "/05_qmmm/19_ensemble_barriers/frame_%s/" % f
        act = [int(x) for x in re.search(r"ActiveAtoms \{([^}]*)\}", open(E + "neb.inp").read()).group(1).split()]
        if act != sorted(set(act)) or len(act) != 1959: stop(rep, "frame %s: neb.inp's ActiveAtoms are not 1959 sorted atoms" % f)
        if act[act.index(QM0):act.index(QM0) + NQM] != list(range(QM0, QM0 + NQM)): stop(rep, "frame %s: QM atoms not all active" % f)
        q = act.index(QM0)
        try:
            M = frames_xyz(dl + "/19_ensemble/frame_%s/neb_MEP.activeRegion_trj.xyz" % f)
            C = frames_xyz(dl + "/19_ensemble/frame_%s/neb_NEB-CI_converged.activeRegion.xyz" % f)
        except (OSError, ValueError) as e:
            stop(rep, "frame %s: %s" % (f, e))
        if len(M) != 10 or len(C) != 1: stop(rep, "frame %s: %d MEP images (10 expected), %d CI structures" % (f, len(M), len(C)))
        want = [SYM.get(z[i], "?") for i in act]
        for c, el, X in M + C:
            if len(el) != 1959: stop(rep, "frame %s: a structure has %d atoms, not 1959" % (f, len(el)))
            if [e.upper() for e in el] != want: stop(rep, "frame %s: elements differ from the prmtop's for the active atoms" % f)
        ref = {"reactant_qm.xyz": M[0][2], "product_qm.xyz": M[-1][2], "transition_state_qm.xyz": C[0][2]}
        dq = {}
        for name, X in ref.items():
            S = frames_xyz(B + name)[0][2]
            dq[name] = float(np.abs(X[q:q + NQM] - S).max())
            if dq[name] > 1e-5: stop(rep, "frame %s: the substrate in the active region differs from %s by %.2e A" % (f, name, dq[name]))
        if f == "41786":
            pc = lambda n: np.array([[float(v) for v in l.split()[1:4]] for l in lzma.open(res + "/phase2.2/cs_stage4d/pc/" + n, "rt").read().splitlines()[1:]])
            row = [i if i < QM0 else i - NQM for i in act if not QM0 <= i < QM0 + NQM]
            mm = [k for k, i in enumerate(act) if not QM0 <= i < QM0 + NQM]
            dR = float(np.abs(M[0][2][mm] - pc("ENZ.pc.xz")[row]).max()); dT = float(np.abs(C[0][2][mm] - pc("ENZTS.pc.xz")[row]).max())
            if dR > 1e-3 or dT > 1e-5: stop(rep, "frame 41786: active MM atoms differ from Stage 4d's env_R by %.2e A, env_TS by %.2e A" % (dR, dT))
            say("frame 41786: the %d active MM atoms equal Stage 4d's env_R (pc/ENZ.pc.xz) within %.1e A and env_TS (pc/ENZTS.pc.xz) within %.1e A"
                % (len(mm), dR, dT))
        d = []
        for nm, seq, cc in RES:
            idx = [act.index(i) for i in range(rp[seq - 1], rp[seq])]
            ci = act.index([i for i in range(rp[seq - 1], rp[seq]) if an[i] == cc][0])
            dd = np.linalg.norm(C[0][2][idx] - M[0][2][idx], axis=1)
            d.append((float(np.linalg.norm(C[0][2][ci] - M[0][2][ci])), float(dd.max())))
        move.append((f, d, max(dq.values())))
    say("all 30 frames: 10 MEP images and the CI, 1959 atoms each, elements as in the prmtop; the 24 substrate atoms of image 0,")
    say("  the last image and the CI equal reactant_qm.xyz, product_qm.xyz and transition_state_qm.xyz (largest |dx| %.1e A)"
        % max(m[2] for m in move))
    say("")
    say("reactant (image 0) -> transition state (CI): movement of each residue, A (charge centre / largest atom)")
    say("frame  " + "  ".join("%-13s" % r[0] for r in RES))
    for f, d, _ in move:
        say(f + "  " + "  ".join("%5.2f / %5.2f  " % x for x in d))
    say("median " + "  ".join("%5.2f / %5.2f  " % (np.median([m[1][k][0] for m in move]), np.median([m[1][k][1] for m in move])) for k in range(len(RES))))
    say("PASS")
    open(rep, "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 5 or sys.argv[1] != "check":
        sys.exit("usage: python3 neb_active_regions.py check <results_repo> <download_folder> <report>")
    check(sys.argv[2], sys.argv[3], sys.argv[4])
