#!/bin/bash
# cs6_d37_harvest.sh - bring DIAGNOSTIC 37's finished runs home, keep them as physical files, read them.
#
#   bash cs6_d37_harvest.sh collect   checks the job has ended, archives every file on hpc1, downloads
#                                     and verifies the archive, places the text outputs in the results
#                                     repository, writes the records and runs the committed reader.
#                                     It commits nothing.
#   bash cs6_d37_harvest.sh commit    commits what collect placed, in both repositories. It does not push.
#
# Rule (Jessica, 18 October 2026): every output file of a finished job is downloaded and kept as a
# physical file. Split decided 28 October 2026 for Diagnostic 14: every file into a sha256-checked
# archive on the PC, the text outputs and a manifest of every archived file into git.
# This script deletes nothing on hpc1. It writes two things there: d37/HPC_MANIFEST.sha256 and the
# archive beside d37.
set -uo pipefail

R="${D37_RESULTS:-$HOME/Desktop/chorismate-thesis-results}"
C="${D37_CODE:-$HOME/Desktop/chorismate-thesis-code}"
A="${D37_ARCHIVE:-$HOME/Desktop/chorismate-hpc-archive}"
H="${D37_HPC:-18660916@hpc}"
U=18660916
S=phase2.2/cs_stage6
W=/home/18660916/system_development/phase2.2/cs_stage6
ARC=hpc_all_files_d37.tar.xz
PIN_INPUTS=9369e1971aebc02cf9920707bfb70140fcd0c4dfb96f7546b16cc10a338bdb5b   # d37/INPUT_CHECKSUMS.sha256
PIN_PREREG=09c83ccee2998cd7851d1f5e4064e7aa80b1ae9379eb5fb750efd05167828627   # STAGE6_DIAGNOSTIC37_20261111.txt
PREREG="$R/$S/STAGE6_DIAGNOSTIC37_20261111.txt"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$HERE/work"
STATE="$HERE/state.env"

say()  { echo "   $*"; }
step() { echo; echo "== $*"; }
die()  { echo; echo "STOPPED: $*"; echo "Nothing has been committed and nothing has been pushed. Nothing on hpc1 was deleted."; exit 1; }
sha()  { sha256sum "$1" | cut -c1-64; }
clean() { [ -z "$(git -C "$1" status --porcelain --untracked-files=all)" ]; }

collect() {
  step "1. This PC"
  for t in git ssh scp tar xz sha256sum python3; do command -v "$t" > /dev/null || die "$t is not installed on this PC."; done
  [ -d "$R/.git" ] || die "$R is not a git repository."
  [ -e "$C/.git" ] || die "$C is not a git repository."
  [ -d "$A" ] || die "the archive folder $A does not exist. Tell Claude where the PC archive lives."
  clean "$R" || die "the results repository has uncommitted or untracked files. If a previous harvest stopped half way, send Claude its output."
  clean "$C" || die "the code repository has uncommitted or untracked files."
  git -C "$R" fetch -q origin || die "could not fetch the results repository."
  git -C "$C" fetch -q origin || die "could not fetch the code repository."
  [ "$(git -C "$R" rev-parse HEAD)" = "$(git -C "$R" rev-parse origin/main)" ] || die "the results repository is not level with GitHub. Pull or push first."
  [ "$(git -C "$C" rev-parse HEAD)" = "$(git -C "$C" rev-parse origin/tier1-realism)" ] || die "the code repository is not level with GitHub. Pull or push first."
  [ -f "$R/$S/d37/INPUT_CHECKSUMS.sha256" ] && [ "$(sha "$R/$S/d37/INPUT_CHECKSUMS.sha256")" = "$PIN_INPUTS" ] || die "the committed d37 inputs are not the ones this harvest was written for."
  [ -f "$PREREG" ] && [ "$(sha "$PREREG")" = "$PIN_PREREG" ] || die "the committed pre-registration is not the one this harvest was written for."
  git -C "$C" cat-file -e HEAD:phase2.2/cs6_d37.py 2> /dev/null || die "the reader cs6_d37.py is not committed in the code repository."
  [ -z "$(find "$R/$S/d37" -name job.out -print -quit)" ] || die "the results repository already holds d37 outputs. This harvest has been run before."
  [ -z "$(ls "$R/$S"/STAGE6_DIAGNOSTIC37_RESULT_*.txt 2> /dev/null)" ] || die "a D37 result document already exists."
  [ ! -e "$A/$S/d37/$ARC" ] || die "$A/$S/d37/$ARC already exists. This harvest has been run before."
  HR="$(git -C "$R" rev-parse HEAD)"; HC="$(git -C "$C" rev-parse HEAD)"
  say "results $HR"; say "code    $HC"; say "archive folder $A"
  rm -rf "$WORK" "$STATE"; mkdir -p "$WORK/x"

  step "2. Has the job ended? (plain qstat output from hpc1 follows)"
  Q="$(ssh "$H" "bash -lc 'qstat -u $U'" 2>&1)" || die "the queue on hpc1 could not be read. A failed query is never taken to mean the job has ended."
  echo "$Q" | sed 's/^/   | /'
  [ -n "$Q" ] || say "| (qstat printed nothing: no jobs of yours are queued or running)"
  echo "$Q" | grep -q "cm_d37" && die "cm_d37 is still in the queue. Run this again when it has gone."
  ssh "$H" "[ -d $W/d37 ] && [ -s $W/d37/run.pbs.out ] && [ ! -e $W/d37/.running ]" \
    || die "on hpc1, d37 is missing, or its run.pbs.out is missing or empty, or its .running lock is still there. The job has not ended cleanly; send Claude the qstat output above."
  ssh "$H" "tail -n 2 $W/d37/run.pbs.out" | sed 's/^/   | /'

  step "3. On hpc1: the inputs that ran, the manifest and the archive (this can take a few minutes)"
  cat > "$WORK/remote.sh" <<EOF
set -euo pipefail
cd $W/d37
sha256sum -c --quiet INPUT_CHECKSUMS.sha256
echo "INPUTS \$(sha256sum INPUT_CHECKSUMS.sha256 | cut -c1-64)"
if [ ! -e HPC_MANIFEST.sha256 ]; then
  find . -type f ! -name HPC_MANIFEST.sha256 ! -name '*.part' | LC_ALL=C sort | sed 's#^\./##' | xargs sha256sum > HPC_MANIFEST.sha256.part
  mv HPC_MANIFEST.sha256.part HPC_MANIFEST.sha256
fi
sha256sum -c --quiet HPC_MANIFEST.sha256
[ "\$(find . -type f ! -name HPC_MANIFEST.sha256 ! -name '*.part' | wc -l)" = "\$(wc -l < HPC_MANIFEST.sha256)" ]
echo "FILES \$(wc -l < HPC_MANIFEST.sha256)"
cd $W
if [ ! -e $ARC ]; then
  tar -cf - d37 | xz -T1 -c > $ARC.part
  mv $ARC.part $ARC
fi
echo "ARCHIVE \$(sha256sum $ARC | cut -c1-64)"
echo "MEMBERS \$(tar -tJf $ARC | grep -vc '/\$')"
EOF
  OUT="$(ssh "$H" bash -s < "$WORK/remote.sh" 2>&1)" || { echo "$OUT" | sed 's/^/   | /'; die "the collection on hpc1 failed. If it names a file that did not match, an input or an output changed on hpc1 after it was written."; }
  RIN="$(echo "$OUT" | awk '$1=="INPUTS"{print $2}')"; NFILES="$(echo "$OUT" | awk '$1=="FILES"{print $2}')"
  RSHA="$(echo "$OUT" | awk '$1=="ARCHIVE"{print $2}')"; NMEM="$(echo "$OUT" | awk '$1=="MEMBERS"{print $2}')"
  [ "$RIN" = "$PIN_INPUTS" ] || die "hpc1 ran inputs whose checksum list is not the committed one."
  [ -n "$NFILES" ] && [ -n "$RSHA" ] && [ -n "$NMEM" ] || die "hpc1's answer could not be read: $OUT"
  [ "$NMEM" = "$((NFILES + 1))" ] || die "the archive on hpc1 holds $NMEM files but the manifest lists $NFILES plus itself."
  say "inputs on hpc1 are the committed ones; $NFILES files in the manifest; archive of $NMEM members"

  step "4. Download, and check the copy against hpc1's"
  mkdir -p "$A/$S/d37" || die "could not create $A/$S/d37"
  scp -q "$H:$W/$ARC" "$A/$S/d37/$ARC.part" || die "the download failed. Nothing was placed."
  [ "$(sha "$A/$S/d37/$ARC.part")" = "$RSHA" ] || { rm -f "$A/$S/d37/$ARC.part"; die "the downloaded archive does not match hpc1's by sha256. The partial copy was removed; run this again."; }
  mv "$A/$S/d37/$ARC.part" "$A/$S/d37/$ARC"
  ABYTES="$(wc -c < "$A/$S/d37/$ARC" | tr -d ' ')"
  say "$A/$S/d37/$ARC  $ABYTES bytes  sha256 $RSHA"

  step "5. Every member of the archive against hpc1's manifest, and its inputs against the repository's"
  tar -xJf "$A/$S/d37/$ARC" -C "$WORK/x" || die "the archive could not be unpacked."
  ( cd "$WORK/x/d37" && sha256sum -c --quiet HPC_MANIFEST.sha256 ) || die "a member of the archive does not match hpc1's manifest."
  [ "$(find "$WORK/x/d37" -type f | wc -l | tr -d ' ')" = "$NMEM" ] || die "the unpacked archive does not hold $NMEM files."
  while read -r _ f; do
    cmp -s "$WORK/x/d37/$f" "$R/$S/d37/$f" || die "$f on hpc1 is not the committed input."
  done < "$R/$S/d37/INPUT_CHECKSUMS.sha256"
  cmp -s "$WORK/x/d37/INPUT_CHECKSUMS.sha256" "$R/$S/d37/INPUT_CHECKSUMS.sha256" || die "INPUT_CHECKSUMS.sha256 on hpc1 is not the committed one."
  say "all $NMEM members match the manifest; all 57 inputs are the committed ones"

  step "6. Place the text outputs and write the records"
  python3 - "$WORK/x/d37" "$R" "$S/d37" "$A" "$S/d37/$ARC" "$ABYTES" "$RSHA" "$NMEM" "$HR" <<'PY' || die "placing the outputs failed; see the line above. Send Claude this output before running anything else."
import hashlib, os, shutil, subprocess, sys
src, repo, rel, arch, arcrel, abytes, asha, nmem, head = sys.argv[1:10]
dst = os.path.join(repo, rel)
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()
committed = set(l.split(None, 1)[1].strip() for l in open(os.path.join(dst, "INPUT_CHECKSUMS.sha256"))) | {"INPUT_CHECKSUMS.sha256"}
files = sorted(os.path.relpath(os.path.join(d, f), src) for d, _, fs in os.walk(src) for f in fs)
inv, new = [], []
for p in files:
    parts = p.split("/")
    if p in committed:
        cat = "in_repo"
    elif (len(parts) == 2 and parts[1] in ("job.out", "job.property.txt")) or p in ("run.pbs.out", "HPC_MANIFEST.sha256"):
        cat = "in_repo"; new.append(p)
    else:
        cat = "pc_only"
    inv.append((os.path.getsize(os.path.join(src, p)), cat, p))
for p in new:
    if os.path.getsize(os.path.join(src, p)) > 95 * 1000 * 1000:
        sys.exit("   REFUSED: %s is over 95 MB and cannot go into git." % p)
    if os.path.exists(os.path.join(dst, p)):
        sys.exit("   REFUSED: %s is already in the repository." % p)
    if subprocess.run(["git", "-C", repo, "check-ignore", "-q", os.path.join(rel, p)]).returncode == 0:
        sys.exit("   REFUSED: %s is covered by .gitignore; git add would skip it silently." % p)
for p in new:
    os.makedirs(os.path.dirname(os.path.join(dst, p)), exist_ok=True)
    shutil.copyfile(os.path.join(src, p), os.path.join(dst, p))
    if sha(os.path.join(src, p)) != sha(os.path.join(dst, p)):
        sys.exit("   REFUSED: %s did not copy byte for byte." % p)
with open(os.path.join(dst, "INVENTORY.tsv"), "w") as f:
    f.write("bytes\tcategory\tpath\n")
    for b, c, p in inv:
        f.write("%d\t%s\t%s\n" % (b, c, p))
inrepo = [p for b, c, p in inv if c == "in_repo"]
with open(os.path.join(dst, "OUTPUTS.sha256"), "w") as f:
    for p in inrepo:
        f.write("%s  %s\n" % (sha(os.path.join(dst, p)), p))
with open(os.path.join(dst, "PC_ARCHIVES.tsv"), "w") as f:
    f.write("# PC archive (relative to the archive folder %s)\tbytes\tsha256\tmembers\n" % arch)
    f.write("%s\t%s\t%s\t%s\n" % (arcrel, abytes, asha, nmem))
runs = sorted(set(p.split("/")[0] for p in committed if "/" in p))
outs = [r for r in runs if os.path.isfile(os.path.join(dst, r, "job.out"))]
norm = [r for r in outs if "ORCA TERMINATED NORMALLY" in open(os.path.join(dst, r, "job.out"), errors="replace").read()]
pc = sorted(set(os.path.basename(p) for b, c, p in inv if c == "pc_only"))
with open(os.path.join(dst, "README.txt"), "w") as f:
    f.write("""d37: the %d runs kept as physical files, then read (cs6_d37_harvest.sh).

Rule (Jessica, 18 October 2026): every output file of a finished job is downloaded and kept as a
physical file, not only reported.

1. PC archive (PC_ARCHIVES.tsv), in the archive folder on Jessica's PC: all %d files of d37 as they
   stood on hpc1 after the job had ended, together with hpc1's own HPC_MANIFEST.sha256, which is
   %s members in one hpc_all_files_d37.tar.xz. hpc1 wrote HPC_MANIFEST.sha256 over its own files
   before the download; the archive was checked against hpc1's copy by sha256 after the download,
   and every archive member was checked against the manifest again where the archive lies.
2. In this repo: job.out and job.property.txt for each run that wrote them, the PBS log
   run.pbs.out, HPC_MANIFEST.sha256, INPUT_CHECKSUMS.sha256, INVENTORY.tsv (every file with its
   size and category) and OUTPUTS.sha256.
3. NOT in this repo: %s.
   Wavefunction and scratch files are excluded by this repository's own .gitignore, the split
   decided on 28 October 2026 for Diagnostic 14: these are fixed-geometry single points whose
   reading uses job.out alone, and the inputs that would regenerate them are committed beside
   them. Every one of these files is in the PC archive with its sha256 in HPC_MANIFEST.sha256.
4. The harvest deletes nothing on hpc1.

job.out is present for %d of the %d runs and reports normal termination in %d. A run that did not
finish is a result and is kept as it is; the reader reports it.

The inputs were committed and pushed before the job was submitted; the harvest started from
results %s. Every input hpc1 ran was checked against the committed one twice: by sha256sum -c on
hpc1 before the archive was made, and byte for byte between the archive's copies and this
repository's.

The archive is named hpc_all_files_d37.tar.xz.
""" % (len(runs), len(files) - 1, nmem, ", ".join(pc) if pc else "nothing", len(outs), len(runs), len(norm), head[:8]))
print("   placed %d files; job.out present for %d of %d runs, normal termination in %d" % (len(new), len(outs), len(runs), len(norm)))
print("   kept in the PC archive only: %s" % (", ".join(pc) if pc else "nothing"))
PY

  step "7. Read: the committed reader, cs6_d37.py"
  RES="STAGE6_DIAGNOSTIC37_RESULT_$(date +%Y%m%d).txt"
  if ! python3 "$C/phase2.2/cs6_d37.py" --res "$R" --code "$C" --prereg "$PREREG" --result "$R/$S/$RES" --tsv "$R/$S/d37.tsv" > "$WORK/reader.out" 2>&1; then
    cat "$WORK/reader.out"
    echo
    echo "THE READER STOPPED AT A GATE, so there is no result document."
    echo "The archive is on your PC and the outputs are in place in the results repository, uncommitted."
    echo "Do NOT run the commit step. Send Claude everything printed above."
    exit 2
  fi
  cat "$WORK/reader.out"
  [ -s "$R/$S/$RES" ] && [ -s "$R/$S/d37.tsv" ] || die "the reader ended without writing its result."

  NMEAS="$(grep -c -E '^   \S.* +measured$' "$R/$S/$RES" || true)"
  NOTM="$(grep -E '^   \S.* +NOT MEASURED$' "$R/$S/$RES" | sed -E 's/^ +//; s/ +NOT MEASURED$//' | paste -sd, - | sed 's/,/, /g' || true)"
  if grep -q '^   NOT EVALUATED' "$R/$S/$RES"; then T2="Test 2, against DIAGNOSTIC 36, is NOT EVALUATED: D36 was not committed when this was read. It is added by a dated addendum."; else T2="Test 2, against DIAGNOSTIC 36, is evaluated in the result document."; fi
  {
    echo "[retained-final] d37: six candidate functionals on the record's own two tests, measured"
    echo
    echo "The control holds: D28's CAM-B3LYP inputs, run again, reproduce D28's energies."
    grep -E '^   def2-(SVPD|TZVPD) +barrier' "$R/$S/$RES" | sed 's/^ *//'
    echo
    echo "$NMEAS of 6 candidates measured.${NOTM:+ Not measured, with the reason in the result document: $NOTM.}"
    echo "$T2"
    echo
    echo "The two tests side by side, as $RES prints them:"
    awk '/^R37.6/{p=1;next} p&&/^   functional/{q=1} q&&NF==0{exit} q{print}' "$R/$S/$RES" | sed 's/^   //'
    echo
    echo "Kept as physical files: job.out and job.property.txt of the runs here, and all $NMEM files in the PC archive (PC_ARCHIVES.tsv, sha256 $RSHA). Nothing was deleted on hpc1. This adopts no level of theory."
    echo
    echo "Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>"
    echo "Claude-Session: https://claude.ai/code/session_01F3vG9e2HqrawvyQuL6J8Ke"
  } > "$HERE/commit_message_results.txt"
  {
    echo "HR=$HR"; echo "HC=$HC"; echo "RES=$RES"; echo "RSHA=$RSHA"; echo "NMEM=$NMEM"
  } > "$STATE"

  step "DONE. Nothing is committed yet."
  say "The result is $R/$S/$RES"
  say "The commit message that the next step will use is $HERE/commit_message_results.txt"
  say "Next: bash $HERE/cs6_d37_harvest.sh commit"
}

commit() {
  [ -s "$STATE" ] && [ -s "$HERE/commit_message_results.txt" ] || die "collect has not finished successfully, so there is nothing to commit."
  # shellcheck disable=SC1090
  . "$STATE"
  [ "$(git -C "$R" rev-parse HEAD)" = "$HR" ] || die "the results repository has moved since collect ran."
  [ "$(git -C "$C" rev-parse HEAD)" = "$HC" ] || die "the code repository has moved since collect ran."
  [ -s "$R/$S/$RES" ] || die "the result document is missing."
  BAD="$(git -C "$R" status --porcelain --untracked-files=all | grep -vc -E "^\?\? $S/(d37/|d37\.tsv$|$RES$)" || true)"
  [ "$BAD" = 0 ] || die "the results repository holds changes that are not part of this harvest."
  [ -n "$(git -C "$R" status --porcelain --untracked-files=all)" ] || die "there is nothing to commit in the results repository."
  clean "$C" || die "the code repository has uncommitted or untracked files."
  ( cd "$R/$S/d37" && sha256sum -c --quiet OUTPUTS.sha256 ) || die "a placed output no longer matches OUTPUTS.sha256."
  cp "$HERE/cs6_d37_harvest.sh" "$C/phase2.2/cs6_d37_harvest.sh" || die "could not copy the harvest script into the code repository."
  git -C "$R" add "$S/d37" "$S/d37.tsv" "$S/$RES" || die "git add failed in the results repository."
  git -C "$R" commit -q -F "$HERE/commit_message_results.txt" || die "the results commit failed."
  git -C "$C" add phase2.2/cs6_d37_harvest.sh || die "git add failed in the code repository."
  git -C "$C" commit -q \
    -m "[retained-final] cs6_d37_harvest.sh: bring DIAGNOSTIC 37 home, keep every file, read it" \
    -m "Checks from plain qstat that the job has ended, has hpc1 write a sha256 manifest of every file and one archive, downloads and verifies the archive against hpc1 and against the manifest, checks the inputs that ran against the committed ones, places job.out and job.property.txt with INVENTORY.tsv, OUTPUTS.sha256, PC_ARCHIVES.tsv and README.txt, and runs the committed cs6_d37.py. It deletes nothing on hpc1 and does not push." \
    -m "Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01F3vG9e2HqrawvyQuL6J8Ke" || die "the code commit failed. The results commit was made; do not push, and send Claude this output."
  clean "$R" && clean "$C" || die "something is left uncommitted after the commits. Do not push; send Claude the output of git status."
  echo
  git -C "$R" log -1 --format='results committed: %h %s'
  git -C "$C" log -1 --format='code committed:    %h %s'
  echo "Nothing has been pushed."
}

case "${1:-}" in
  collect) collect ;;
  commit)  commit ;;
  *) echo "usage: bash cs6_d37_harvest.sh collect | commit"; exit 1 ;;
esac
