#!/usr/bin/env python3
"""
cs6_prepare.py - Stage 6 of the charged-sphere development (phase2.2, cs_stage6/STAGE6_CRITERIA.txt with
STAGE6_AMENDMENT1.txt): the substrate geometries Stage 6 relaxes, taken from frame 42436's committed smoke-test
stationary points.

  python3 cs6_prepare.py extract <results_repo> <OUT>

EXTRACT
Reads, from <results_repo>/05_qmmm/24_sp_barriers/f42436_r8_a1/, the three committed archives
hpc_text_outputs_2_{ts,r_irc,p_irc}.tar.xz, and from each the step's final.xyz and STAGE_RESULT.txt. Those archives
are what SP_AMENDMENT2.txt records as frame 42436's package run (smoke test 2), and STAGE5_CRITERIA.txt section 8
names the same three final.xyz as Stage 5 production's per-frame inputs.

Checks, all of them from STAGE6_CRITERIA.txt and cs5_align.py's extract, and every one a stop:
  (a) each archive's sha256 equals the value STAGE6_CRITERIA.txt records;
  (b) each step's STAGE_RESULT.txt has status=PASS;
  (c) each final.xyz's sha256 equals that record's xyz_sha256;
  (d) all three full-system files have the same atom count;
  (e) the 24 substrate atoms are ORCA's QMAtoms {6207:6230} (0-based) and their elements are the substrate's own
      order, C H H C C O O O C H C H C C H C H C H O H C O O;
  (f) the breaking O3-C4 and forming C1-C6 distances of the extracted substrate equal the record's break_A and
      form_A to 3 decimals;
  (g) each step's recorded "substrate" character is the one that step should have: reactant for r_irc, neither for
      ts, product for p_irc.

Writes <OUT>/frame_42436/{reactant_qm.xyz, transition_state_qm.xyz, product_qm.xyz} - the names cs5_grid.py and
cs5_amatrix.py read - plus <OUT>/frame_42436/transform.txt, which cs5_amatrix.py requires and which here is the
IDENTITY (rotation I, translation 0), in cs5_align.py's own layout and read by its own np.loadtxt call: rows 0-2 the
rotation, row 3 the translation. The coordinates are copied from the source files verbatim, at their full 14 decimal
places, rather than reformatted, so the extracted geometry is byte-traceable to the committed final.xyz. And
<OUT>/EXTRACT.txt recording every source, sha256 and checked value. Refuses to write into
an existing folder. Deterministic: a second run into another folder is byte-identical. Writes nothing else and runs
nothing.

NO ALIGNMENT IS APPLIED, and that is not an omission. Stage 5 fits one rotation and translation per frame on the
rigid core of that frame's R onto the REFERENCE frame's R and applies it unchanged to that frame's TS and P
(cs5_align.py align), so that the reactant-to-TS motion is preserved. With a single frame there is no other frame to
align to, and a rigid motion applied to all three states together cannot change the grid relative to the substrate -
only the coordinate values, which the in vacuo energy and potential do not depend on (cs5_amatrix.py's check (b)
measures exactly that). The three states are already in one coordinate frame, the full system's, because all three
optimisations ran in it.

Python 3.6, standard library only.
"""
import hashlib, lzma, math, os, re, sys, tarfile

NAMES = "C1 H1 H2 C2 C3 O1 O2 O3 C4 H3 C5 H4 C6 C7 H5 C8 H6 C9 H7 O4 H8 C10 O5 O6".split()
ELS = [re.sub(r"\d+$", "", n) for n in NAMES]
QM0, NQM = 6207, 24
BREAK = (NAMES.index("O3"), NAMES.index("C4"))
FORM = (NAMES.index("C1"), NAMES.index("C6"))
FRAME = "42436"
RUN = "f42436_r8_a1"
STATES = [("r_irc", "reactant_qm.xyz", "reactant"), ("ts", "transition_state_qm.xyz", "neither"),
          ("p_irc", "product_qm.xyz", "product")]
# sha256 of the three archives, as STAGE6_CRITERIA.txt records them
ARCHIVES = {
    "ts": "6f01f6a3efe5c0d0b94c0de9de805a23cc126e7fde7e9f50de25bdbf465c9022",
    "r_irc": "3f9ae488565c11d6b9aa14f98886a5e37e6dfe0b4e314038778a35d359bf354b",
    "p_irc": "5a7dc261ad844bcb80e1679a120595ae84a7bfae5b5be1f9bec3f2ba6817d326",
}

def stop(m):
    sys.stderr.write("STOP: " + m + "\n")
    sys.exit(1)

def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()

def sha_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

def member(tf, suffix):
    hits = [m for m in tf.getmembers() if m.isfile() and m.name.endswith(suffix)]
    if len(hits) != 1:
        stop("%s: expected exactly one member ending %s, found %d" % (tf.name, suffix, len(hits)))
    return tf.extractfile(hits[0]).read()

def record(text):
    out = {}
    for ln in text.splitlines():
        if "=" in ln:
            k, v = ln.split("=", 1)
            out[k.strip()] = v.strip()
    return out

def extract(repo, out):
    src = os.path.join(repo, "05_qmmm", "24_sp_barriers", RUN)
    if not os.path.isdir(src):
        stop("%s is not there: is this the results repo?" % src)
    if os.path.exists(out):
        stop("%s exists already; this script never writes into an existing folder" % out)
    log = ["cs6_prepare.py extract - Stage 6 substrate geometries (STAGE6_CRITERIA.txt, STAGE6_AMENDMENT1.txt)",
           "source: %s" % os.path.join("05_qmmm", "24_sp_barriers", RUN), ""]
    got = {}
    counts = set()
    for step, outname, want_char in STATES:
        arch = os.path.join(src, "hpc_text_outputs_2_%s.tar.xz" % step)
        if not os.path.isfile(arch):
            stop("%s is missing" % arch)
        h = sha_file(arch)
        if h != ARCHIVES[step]:
            stop("%s: sha256 %s, but STAGE6_CRITERIA.txt records %s" % (os.path.basename(arch), h, ARCHIVES[step]))
        with tarfile.open(arch, "r:xz") as tf:
            res_b = member(tf, "/%s/STAGE_RESULT.txt" % step)
            xyz_b = member(tf, "/%s/final.xyz" % step)
        r = record(res_b.decode())
        if r.get("status") != "PASS":
            stop("%s: STAGE_RESULT.txt status is %r, not PASS" % (step, r.get("status")))
        hx = sha_bytes(xyz_b)
        if hx != r.get("xyz_sha256"):
            stop("%s: final.xyz sha256 %s, but its STAGE_RESULT.txt records %s" % (step, hx, r.get("xyz_sha256")))
        lines = xyz_b.decode().splitlines()
        n = int(lines[0].split()[0])
        counts.add(n)
        body = lines[2:2 + n]
        if len(body) != n:
            stop("%s: final.xyz says %d atoms but holds %d lines" % (step, n, len(body)))
        sub = body[QM0:QM0 + NQM]
        if len(sub) != NQM:
            stop("%s: fewer than %d atoms from index %d" % (step, NQM, QM0))
        els = [l.split()[0] for l in sub]
        if els != ELS:
            stop("%s: the 24 atoms at %d are %s, not the substrate's order %s"
                 % (step, QM0, " ".join(els), " ".join(ELS)))
        xyz = [tuple(float(x) for x in l.split()[1:4]) for l in sub]
        def dist(a, b):
            return math.sqrt(sum((xyz[a][k] - xyz[b][k]) ** 2 for k in range(3)))
        br, fo = dist(*BREAK), dist(*FORM)
        for label, value, key in (("O3-C4", br, "break_A"), ("C1-C6", fo, "form_A")):
            want = float(r[key])
            if abs(round(value, 3) - want) > 1e-9:
                stop("%s: %s is %.3f but STAGE_RESULT.txt records %s=%s" % (step, label, value, key, r[key]))
        if r.get("substrate") != want_char:
            stop("%s: STAGE_RESULT.txt calls the substrate %r, expected %r" % (step, r.get("substrate"), want_char))
        text = "%d\n%s: frame %s %s, the 24 substrate atoms of %s/%s/final.xyz\n%s\n" % (
            NQM, "cs6_prepare.py", FRAME, step, RUN, step, "\n".join(sub))
        got[outname] = text
        log.append("%-22s from %s/final.xyz" % (outname, step))
        log.append("    archive sha256   %s" % ARCHIVES[step])
        log.append("    final.xyz sha256 %s  (= STAGE_RESULT.txt xyz_sha256, status PASS, %s cycles)"
                   % (hx, r.get("cycles")))
        log.append("    atoms %d in the full system; substrate = indices %d-%d, elements checked"
                   % (n, QM0, QM0 + NQM - 1))
        log.append("    O3-C4 %.3f A and C1-C6 %.3f A, equal to break_A and form_A; substrate=%s"
                   % (br, fo, r.get("substrate")))
        log.append("    energy_Eh %s" % r.get("energy_Eh"))
        log.append("")
    if len(counts) != 1:
        stop("the three full-system files have different atom counts: %s" % sorted(counts))
    log.append("no transform is applied: one frame, so there is no reference to align to, and a rigid motion of all")
    log.append("three states together cannot change the grid relative to the substrate (see this script's docstring)")
    log.append("")
    d = os.path.join(out, "frame_" + FRAME)
    os.makedirs(d)
    # cs5_amatrix.py reads frame_<F>/transform.txt (np.loadtxt: rows 0-2 rotation, row 3 translation). With one
    # frame no alignment is applied, so it is the identity; writing it explicitly records that, and keeps
    # cs5_amatrix.py's inverse-transform path exact (original = (aligned - 0) @ I.T = aligned).
    tf = ["# rows 0-2: rotation R; row 3: translation t",
          "# aligned = original @ R + t ; original = (aligned - t) @ R.T",
          "# the IDENTITY: cs6_prepare.py applies no alignment, because Stage 6 uses one frame and a rigid motion",
          "# of all three states together cannot change the grid relative to the substrate"]
    for row in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0)):
        tf.append(" ".join("%.12f" % v for v in row))
    with open(os.path.join(d, "transform.txt"), "w") as f:
        f.write("\n".join(tf) + "\n")
    log.append("wrote frame_%s/transform.txt        the identity, as cs5_amatrix.py requires" % FRAME)
    for name, text in sorted(got.items()):
        with open(os.path.join(d, name), "w") as f:
            f.write(text)
        log.append("wrote frame_%s/%-22s sha256 %s" % (FRAME, name, sha_bytes(text.encode())))
    with open(os.path.join(out, "EXTRACT.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print("\n".join(log))

def main():
    if len(sys.argv) != 4 or sys.argv[1] != "extract":
        sys.stderr.write(__doc__)
        sys.exit(2)
    extract(sys.argv[2], sys.argv[3])

main()
