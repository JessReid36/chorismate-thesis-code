#!/usr/bin/env python3
"""cs6_d39_prepare.py  DIAGNOSTIC 39: the detachment energy at the benchmark level.

Generates six single-point inputs from Diagnostic 38's committed d38/N/job.inp by
applying only the substitutions the pre-registration documents, and verifies each one by
rebuilding it from the template rather than by trusting its own substitution.

Nothing about the charge arrangements is re-derived here. D38 derived them from each
design's committed site.ORCAFF.prms and combined.xyz and its gates checked them; this
script copies d38/B1/charge.pc and d38/A1/charge.pc byte for byte and refuses if a copy
differs from its original.

Usage: cs6_d39_prepare.py <results-repo> <out-dir> <hpc-dir>
"""
import hashlib, io, os, shutil, sys

S6 = "phase2.2/cs_stage6"
HPCUSER = "18660916"
MARK = "d39"

# run name -> (route, source d38 run whose charge file and dianion output it pairs with)
RUNS = (
    ("M_N",  "M", None),
    ("M_B1", "M", "B1"),
    ("M_A1", "M", "A1"),
    ("X_N",  "X", None),
    ("X_B1", "X", "B1"),
    ("X_A1", "X", "A1"),
)

# The three dianion halves D39 reuses instead of re-running. sha256 of the committed
# outputs, recorded in section 6 of the pre-registration.
DIANION_SHA = {
    "N":  "18a58ac48f7890b6bbd7622a507cd560a75ce8f13cb3742fa6cec9ce6ef28a9c",
    "B1": "307eab716f3f719ebff3afc9d07b58dbfc891abb7ace04599b6fcb04e6835929",
    "A1": "64013339b8e8f5ee8442638f785cebd32434512cec7a4551335ce796d8228422",
}

# scratch floor, set from the corrected figure in
# STAGE6_DIAGNOSTIC38_RESULT_ERRATUM1_20261115.txt and not from the 11 MB the result
# reported. 31 GB <= closed-shell peak < 43 GB; the manual puts open-shell at about
# twice the cost; this is above twice the corrected upper bound.
WANT_GB = 150
MIN_GB = 100
SAMPLE_S = 60


def die(msg):
    sys.stderr.write("REFUSED: %s\n" % msg)
    raise SystemExit(1)


def rd(p):
    return io.open(p, encoding="utf-8").read()


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()


def build_input(tmpl_lines, route, charged):
    """Apply only the documented substitutions to D38's template.

    route M : the charge/multiplicity line becomes -1 2.  Nothing else changes.
    route X : the method keyword becomes IP-EOM-DLPNO-CCSD and an %mdci nroots block
              is inserted after %scf.  The charge/multiplicity line is untouched.
    charged : a %pointcharges line is inserted immediately before the geometry line.
    """
    out = []
    for ln in tmpl_lines:
        if route == "X" and ln.startswith("!"):
            if "DLPNO-CCSD(T)" not in ln:
                die("the template's keyword line does not contain DLPNO-CCSD(T): %r" % ln)
            ln = ln.replace("DLPNO-CCSD(T)", "IP-EOM-DLPNO-CCSD", 1)
        if ln.startswith("*"):
            if charged:
                out.append('%pointcharges "charge.pc"')
            if route == "M":
                if ln != "* xyzfile -2 1 substrate.xyz":
                    die("the template's geometry line is not the expected dianion line: %r" % ln)
                ln = "* xyzfile -1 2 substrate.xyz"
        out.append(ln)
        if route == "X" and ln.startswith("%scf"):
            out.append("%mdci nroots 2 end")
    return "\n".join(out) + "\n"


def jobscript(hpc, runs):
    nm = " ".join(r for r, _, _ in runs)
    L = []
    a = L.append
    a("#!/bin/bash")
    a("#PBS -N cm_d39")
    a("#PBS -l select=1:ncpus=16:mem=120gb")
    a("#PBS -l walltime=24:00:00")
    a("#PBS -m ae")
    a("#PBS -M %s@sun.ac.za" % HPCUSER)
    a("#PBS -j oe")
    a("#PBS -o %s/run.pbs.out" % hpc)
    a("# DIAGNOSTIC 39: six single points on one geometry. Three are monoanions at")
    a("#   charge -1 multiplicity 2 for a CCSD(T) detachment energy against D38's")
    a("#   committed dianions; three compute the detachment directly from the closed-")
    a("#   shell dianion with IP-EOM-DLPNO-CCSD, as a control on the open-shell")
    a("#   reference the first three use.")
    a("# Sequential, one run's scratch cleared before the next starts.")
    a("# walltime 24:00:00 routes to day automatically. No -q.")
    a("#")
    a("# SCRATCH. D38's result reported its peak as 11 MB and said to size later jobs")
    a("#   from it. That measurement was taken after ORCA returned, so it recorded the")
    a("#   residue. d38/PROGRESS.txt, same commit, records used=31G. The erratum states")
    a("#   the corrected bound, 31 GB <= peak < 43 GB closed-shell, and this script is")
    a("#   sized from that: %d GB asked for, %d GB the least a candidate may offer."
      % (WANT_GB, MIN_GB))
    a("#   That rules out $TMPDIR here (38-43 GB) on purpose. It also samples the peak")
    a("#   every %d s WHILE ORCA runs, instead of once after it has finished." % SAMPLE_S)
    a('export PATH="/apps/openmpi/4.1.1/bin:$PATH"')
    a('export LD_LIBRARY_PATH="/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:'
      '/apps/mambaforge/envs/medaka/lib:${LD_LIBRARY_PATH:-}"')
    a("export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1")
    a("ORCA=/home/apps2/ORCA/6.0.1")
    a("H=%s" % hpc)
    a('P="$H/PROGRESS.txt"')
    a('PK="$H/.peak_mb"')
    a('RUNS="%s"' % nm)
    a("")
    a('cd "$H" || exit 1')
    a('if mkdir .running 2>/dev/null; then :; else')
    a('  echo "FAIL: another job holds .running in $H"; exit 1')
    a("fi")
    a("")
    a('MARK="%s.${PBS_JOBID:-nojobid}"' % MARK)
    a("WANT_GB=%d" % WANT_GB)
    a("MIN_GB=%d" % MIN_GB)
    a('S=""')
    a('AVAIL=""')
    a('SMARGIN=""')
    a("RESIDUE=0")
    a('echo 0 > "$PK"')
    a("")
    a("pick() {")
    a('  for c in /scratch-small-local/%s /scratch/%s "${TMPDIR:-}"; do' % (HPCUSER, HPCUSER))
    a('    if [ -z "$c" ]; then continue; fi')
    a('    t="$c/$MARK"')
    a('    if mkdir -p "$t" 2>/dev/null; then :; else')
    a('      echo "   candidate $c  mkdir failed"; continue; fi')
    a('    if echo ok > "$t/.probe" 2>/dev/null; then :; else')
    a('      echo "   candidate $c  created but not writable"; continue; fi')
    a('    rm -f "$t/.probe"')
    a("    a=$(df -P -BG \"$t\" 2>/dev/null | awk 'NR==2{sub(/G$/,\"\",$4); print $4}')")
    a('    case "${a:-}" in')
    a('      \'\'|*[!0-9]*) echo "   candidate $c  free space unreadable"; continue ;;')
    a("    esac")
    a('    if [ "$a" -ge "$1" ]; then S="$t"; AVAIL="$a"; return 0; fi')
    a('    echo "   candidate $c  ${a}G free, under the ${1}G asked for"')
    a("  done")
    a("  return 1")
    a("}")
    a("")
    a("cleanup() {")
    a('  kill "${WATCH:-}" 2>/dev/null')
    a('  if [ -n "$S" ] && [ -d "$S" ]; then')
    a("    for r in $RUNS; do")
    a('      if [ -f "$S/$r/job.out" ] && [ -d "$H/$r" ]; then')
    a('        cp -f "$S/$r/job.out" "$H/$r/job.out" 2>/dev/null')
    a("      fi")
    a("    done")
    a('    echo "scratch at exit: $(du -sh "$S" 2>/dev/null | cut -f1)"')
    a('    case "$S" in')
    a('      */%s.*) rm -rf "$S" ;;' % MARK)
    a('      *) echo "REFUSING to remove an unexpected scratch path: $S" ;;')
    a("    esac")
    a("  fi")
    a('  rmdir "$H/.running" 2>/dev/null')
    a("}")
    a("trap cleanup EXIT")
    a("")
    a('echo "host=$(hostname) start=$(date) jobid=${PBS_JOBID:-none}"')
    a('echo "cores=$(nproc) ; memory:"; free -g | head -2')
    a("")
    a('echo "choosing scratch, in order of preference:"')
    a('if pick "$WANT_GB"; then')
    a('  SMARGIN="full"')
    a('elif pick "$MIN_GB"; then')
    a('  SMARGIN="REDUCED"')
    a('  echo "   NOTE: reduced margin. The corrected closed-shell figure is 31 GB to"')
    a('  echo "   43 GB per run and the manual puts open-shell at about twice that, so"')
    a('  echo "   this is near the edge. R39.9 reports the sampled maximum."')
    a("else")
    a('  echo "FAIL: no candidate on $(hostname) is writable with at least ${MIN_GB}G free."')
    a('  echo "      This floor comes from a measurement, not from caution: D38 needed at"')
    a('  echo "      least 31 GB per closed-shell run and these runs are open-shell."')
    a('  echo "      Nothing is wrong with the calculation. Resubmit and PBS will usually"')
    a('  echo "      place the job on another node. $TMPDIR is deliberately too small to"')
    a('  echo "      be accepted here, so a node whose /scratch-small-local is full will"')
    a('  echo "      always refuse."')
    a("  exit 1")
    a("fi")
    a('echo "R39.9 scratch=$S avail=${AVAIL}G margin=$SMARGIN"')
    a("")
    a("printf '%s  d39 on %s, scratch %s, %sG free, sampling every %ds\\n' \\")
    a("  \"$(date '+%F %T')\" \"$(hostname)\" \"$S\" \"$AVAIL\" " + str(SAMPLE_S) + " > \"$P\"")
    a("")
    a("# R39.9. The watcher carries the maximum it sees while ORCA is running and writes")
    a("# it to $PK, because a measurement taken after ORCA returns records the residue.")
    a("# That was Diagnostic 38's defect and its erratum records it.")
    a("( sleep 20")
    a("  M=0")
    a('  while [ -d "$S" ]; do')
    a('    u=$(du -sm "$S" 2>/dev/null | cut -f1)')
    a('    case "${u:-}" in')
    a("      ''|*[!0-9]*) u=0 ;;")
    a("    esac")
    a('    if [ "$u" -gt "$M" ]; then M="$u"; echo "$M" > "$PK"; fi')
    a("    printf '%s  used=%sM  max=%sM  storage=%s  running=%s\\n' \\")
    a('      "$(date \'+%F %T\')" "$u" "$M" \\')
    a('      "$(grep -hm1 \'PNO-Integral Storage\' "$S"/*/job.out 2>/dev/null '
      "| sed 's/.*\\.\\.\\.[[:space:]]*//' | head -1)\" \\")
    a('      "$(ls "$S" 2>/dev/null | tr \'\\n\' \' \')" >> "$P"')
    a("    sleep %d" % SAMPLE_S)
    a("  done ) &")
    a("WATCH=$!")
    a("")
    a("RC=0")
    a("SKIP_X=0")
    a("for r in $RUNS; do")
    a('  if [ -s "$H/$r/job.out" ] && grep -q "ORCA TERMINATED NORMALLY" "$H/$r/job.out"; then')
    a('    echo "=== $r done before, skipping"; continue')
    a("  fi")
    a("  # C39.2. If X_N does not run, the other two X runs are not attempted and the")
    a("  # result records the refusal. Route M then stands alone, uncontrolled.")
    a('  case "$r" in')
    a("    X_*)")
    a('      if [ "$SKIP_X" -eq 1 ]; then')
    a('        echo "=== $r NOT ATTEMPTED: X_N did not terminate normally (C39.2)"')
    a("        printf '%s  %s not attempted, C39.2\\n' \"$(date '+%F %T')\" \"$r\" >> \"$P\"")
    a("        continue")
    a("      fi")
    a("      ;;")
    a("  esac")
    a('  echo "=== $r  start $(date)"')
    a('  mkdir -p "$S/$r" || { RC=1; break; }')
    a('  cp -p "$H/$r/job.inp" "$H/$r/substrate.xyz" "$S/$r/" || { RC=1; break; }')
    a('  if [ -f "$H/$r/charge.pc" ]; then cp -p "$H/$r/charge.pc" "$S/$r/"; fi')
    a('  ( cd "$S/$r" && "$ORCA/orca" job.inp > job.out 2>&1 </dev/null )')
    a("  rc=$?")
    a('  cp -f "$S/$r/job.out" "$H/$r/job.out" 2>/dev/null')
    a("  for f in job.property.txt job.gbw; do")
    a('    if [ -f "$S/$r/$f" ]; then cp -f "$S/$r/$f" "$H/$r/$f"; fi')
    a("  done")
    a('  echo "$r orca_exit=$rc normal=$(grep -c \'ORCA TERMINATED NORMALLY\' "$H/$r/job.out")"')
    a('  echo "$r storage: $(grep -m1 \'PNO-Integral Storage\' "$H/$r/job.out")"')
    a('  echo "$r E(0):    $(grep -m1 \'E(0)\' "$H/$r/job.out")"')
    a('  echo "$r E(CCSD): $(grep -m1 \'E(CCSD)\' "$H/$r/job.out")"')
    a('  echo "$r CCSD(T): $(grep -m1 \'E(CCSD(T))\' "$H/$r/job.out")"')
    a('  echo "$r S**2:    $(grep -m1 \'Expectation value of <S\\*\\*2>\' "$H/$r/job.out")"')
    a('  echo "$r IProot:  $(grep -m1 \'IROOT\' "$H/$r/job.out")"')
    a('  grep -m1 "TOTAL RUN TIME" "$H/$r/job.out"')
    a("  printf '%s  %s orca_exit=%s\\n' \"$(date '+%F %T')\" \"$r\" \"$rc\" >> \"$P\"")
    a("  # C39.2's trigger, read from the output and not from the exit code alone.")
    a('  if [ "$r" = "X_N" ]; then')
    a('    if grep -q "ORCA TERMINATED NORMALLY" "$H/$r/job.out"; then :; else')
    a("      SKIP_X=1")
    a('      echo "X_N did not terminate normally. X_B1 and X_A1 will not be attempted."')
    a('      echo "The reason, as ORCA states it:"')
    a('      grep -iE "aborting|error|not available|unknown|no auxiliary|basis" '
      '"$H/$r/job.out" | tail -20')
    a("    fi")
    a("  fi")
    a("  # The residue, kept under its own name. This is NOT the peak; $PK holds that.")
    a('  u=$(du -sm "$S" 2>/dev/null | cut -f1)')
    a('  case "${u:-}" in')
    a("    ''|*[!0-9]*) u=0 ;;")
    a("  esac")
    a('  if [ "$u" -gt "$RESIDUE" ]; then RESIDUE="$u"; fi')
    a('  echo "$r residue after this run: ${u} MB  (not the peak; see R39.9)"')
    a('  case "$S" in')
    a('    */%s.*) rm -rf "$S/$r" ;;' % MARK)
    a('    *) echo "REFUSING to remove an unexpected scratch path: $S/$r" ;;')
    a("  esac")
    a('  if [ "$rc" -ne 0 ]; then RC=$rc; fi')
    a("done")
    a('PEAK=$(cat "$PK" 2>/dev/null || echo 0)')
    a('echo "R39.9 sampled maximum while running: ${PEAK} MB, every %ds, on $S at the '
      '$SMARGIN margin"' % SAMPLE_S)
    a('echo "R39.9 largest residue after a run:  ${RESIDUE} MB"')
    a('echo "R39.9 the sampled maximum is a lower bound on the true peak"')
    a('echo "all runs finished, worst exit=$RC $(date)"')
    a("exit $RC")
    return "\n".join(L) + "\n"


def main():
    if len(sys.argv) != 4:
        sys.stderr.write("usage: cs6_d39_prepare.py <results-repo> <out-dir> <hpc-dir>\n")
        return 2
    res, out, hpc = sys.argv[1:4]
    s6 = os.path.join(res, S6)
    d38 = os.path.join(s6, "d38")

    tp = os.path.join(d38, "N", "job.inp")
    if not os.path.isfile(tp):
        die("%s is missing. Every D39 input is built from it." % tp)
    tmpl = rd(tp)
    tl = tmpl.rstrip("\n").split("\n")
    if len(tl) != 5:
        die("the template is %d lines, expected 5: %r" % (len(tl), tl))
    for want, got in (("!", tl[0]), ("%maxcore", tl[1]), ("%pal", tl[2]),
                      ("%scf", tl[3]), ("*", tl[4])):
        if not got.startswith(want):
            die("template line does not start with %r: %r" % (want, got))
    print("template  %s" % tp)
    for ln in tl:
        print("          %s" % ln)

    geo = os.path.join(s6, "d33a", "sp_42436_R", "substrate.xyz")
    if not os.path.isfile(geo):
        die("%s is missing. It is the geometry every run uses." % geo)
    gsha = sha(geo)
    if sha(os.path.join(d38, "N", "substrate.xyz")) != gsha:
        die("d38/N/substrate.xyz is not the committed d33a geometry. Refusing.")
    print("geometry  %s  sha256 %s" % (geo, gsha[:16]))

    for k, want in sorted(DIANION_SHA.items()):
        p = os.path.join(d38, k, "job.out")
        if not os.path.isfile(p):
            die("%s is missing. D39 reuses it instead of re-running the dianion." % p)
        got = sha(p)
        if got != want:
            die("d38/%s/job.out has sha256 %s, section 6 records %s.\n"
                "   The dianion half of a subtraction has changed. Refusing." % (k, got, want))
    print("dianions  d38/N, d38/B1, d38/A1 job.out all match their recorded sha256")

    if os.path.exists(out):
        shutil.rmtree(out)
    os.makedirs(out)

    manifest = []
    print()
    for name, route, src in RUNS:
        rd_ = os.path.join(out, name)
        os.makedirs(rd_)
        charged = src is not None
        txt = build_input(tl, route, charged)
        # Verify by reconstruction: rebuild from the template and require equality.
        if build_input(rd(tp).rstrip("\n").split("\n"), route, charged) != txt:
            die("%s does not reproduce from the template. Refusing." % name)
        io.open(os.path.join(rd_, "job.inp"), "w", encoding="utf-8").write(txt)
        shutil.copyfile(geo, os.path.join(rd_, "substrate.xyz"))
        if sha(os.path.join(rd_, "substrate.xyz")) != gsha:
            die("%s/substrate.xyz did not copy byte for byte." % name)
        if charged:
            cp = os.path.join(d38, src, "charge.pc")
            if not os.path.isfile(cp):
                die("%s is missing." % cp)
            shutil.copyfile(cp, os.path.join(rd_, "charge.pc"))
            if sha(os.path.join(rd_, "charge.pc")) != sha(cp):
                die("%s/charge.pc did not copy byte for byte from d38/%s." % (name, src))
        # Every line that is not the method line, the geometry line, the inserted
        # %mdci line or the inserted %pointcharges line must be the template's, char
        # for char.
        got = txt.rstrip("\n").split("\n")
        keep_t = [l for l in tl if not (l.startswith("!") or l.startswith("*"))]
        keep_g = [l for l in got if not (l.startswith("!") or l.startswith("*")
                                         or l.startswith("%mdci")
                                         or l.startswith("%pointcharges"))]
        if keep_g != keep_t:
            die("%s changed a line it must not: %r vs %r" % (name, keep_g, keep_t))
        gl = [l for l in got if l.startswith("*")]
        if len(gl) != 1:
            die("%s has %d geometry lines." % (name, len(gl)))
        want_geo = ("* xyzfile -1 2 substrate.xyz" if route == "M"
                    else "* xyzfile -2 1 substrate.xyz")
        if gl[0] != want_geo:
            die("%s geometry line is %r, expected %r" % (name, gl[0], want_geo))
        npc = 1 if charged else 0
        if len([l for l in got if l.startswith("%pointcharges")]) != npc:
            die("%s has the wrong number of %%pointcharges lines." % name)
        if route == "X":
            if "IP-EOM-DLPNO-CCSD" not in got[0]:
                die("%s is a route X run without the IP-EOM keyword." % name)
            if "%mdci nroots 2 end" not in got:
                die("%s is a route X run without the nroots block." % name)
        else:
            if got[0] != tl[0]:
                die("%s is a route M run whose keyword line is not the template's." % name)
            if any(l.startswith("%mdci") for l in got):
                die("%s is a route M run carrying an %%mdci block." % name)
        for f in ("job.inp", "substrate.xyz") + (("charge.pc",) if charged else ()):
            manifest.append(("%s/%s" % (name, f), sha(os.path.join(rd_, f))))
        print("  %-5s route %s  charge/mult %s  %s"
              % (name, route, want_geo.split()[2] + " " + want_geo.split()[3],
                 ("charge.pc from d38/%s" % src) if charged else "no embedding charge"))

    js = jobscript(hpc, RUNS)
    for bad in ("#PBS -q", "qsub "):
        if bad in js:
            die("the generated job script contains %r. Allocation is automatic." % bad)
    if "du -sm" not in js or "$PK" not in js:
        die("the generated job script has no live peak sampler.")
    io.open(os.path.join(out, "run.pbs"), "w", encoding="utf-8").write(js)

    mp = os.path.join(out, "INPUT_CHECKSUMS.sha256")
    with io.open(mp, "w", encoding="utf-8") as f:
        for n, h in manifest:
            f.write("%s  %s\n" % (h, n))
    print()
    print("  run.pbs                 16 cores, 120 gb, 24:00:00, no queue")
    print("  INPUT_CHECKSUMS.sha256  %d inputs" % len(manifest))
    print("  scratch floor           %d GB asked for, %d GB minimum, sampled every %d s"
          % (WANT_GB, MIN_GB, SAMPLE_S))
    return 0


if __name__ == "__main__":
    sys.exit(main())
