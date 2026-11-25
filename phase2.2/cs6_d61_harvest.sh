#!/bin/bash
# cs6_d61_harvest.sh - bring the finished runs of STAGE 6, DIAGNOSTIC 61 home, keep them as physical
# files, and read them.
#
#   bash cs6_d61_harvest.sh collect
#        checks that every job has ended and that its log reaches its last line, archives
#        every file on hpc1, downloads and verifies the archive, places the text outputs in the
#        results repository, writes the records and runs the committed reader. It commits nothing.
#   bash cs6_d61_harvest.sh collect anyway
#        the same, for jobs that were cut short or never submitted and that are NOT going to be
#        submitted again. The runs that never ran are then reported NOT RUN by the reader.
#   bash cs6_d61_harvest.sh commit
#        commits what collect placed, in the results repository. It does not push.
#
# The committed cs6_d59_harvest.sh, with what names the diagnostic changed: DIAGNOSTIC 61 is 13
# jobs (J01 to J13), each with its job script run_<job>.pbs, its log run_<job>.pbs.out and its lock
# .running_<job>; the text files kept in the repository are job.out and job.property.txt; the
# reader writes two tables; and step 1 also checks that D59's result is committed, as D61's reader
# reads D59's runs N and GN. No job script makes a run a second time in its folder (D44
# NOTE 3). The archive folder is on Jessica's external drive and is given as CS6_ARCHIVE; the work
# folder, which holds the unpacked check copy, can be given as CS6_WORK.
# Rule (Jessica, 18 October 2026): every output file of a finished job is downloaded and kept as a
# physical file. Split decided 28 October 2026 for Diagnostic 14: every file into a sha256-checked
# archive on the PC, the text outputs and a manifest of every archived file into git.
# This script deletes nothing, on hpc1 or on the PC. On hpc1 it writes two things:
# d61/HPC_MANIFEST.sha256 and the archive beside d61.
set -uo pipefail

R="${CS6_RESULTS:-$HOME/Desktop/chorismate-thesis-results}"
C="${CS6_CODE:-$HOME/Desktop/chorismate-thesis-code}"
A="${CS6_ARCHIVE:-$HOME/Desktop/chorismate-hpc-archive}"
H="${CS6_HPC:-18660916@hpc}"
U=18660916
S=phase2.2/cs_stage6
W=/home/18660916/system_development/phase2.2/cs_stage6

DIR=d61; NIN=196; NRUN=120
LEVELS="J01 J02 J03 J04 J05 J06 J07 J08 J09 J10 J11 J12 J13"
JOBS="cm_d61_j01|cm_d61_j02|cm_d61_j03|cm_d61_j04|cm_d61_j05|cm_d61_j06|cm_d61_j07|cm_d61_j08|cm_d61_j09|cm_d61_j10|cm_d61_j11|cm_d61_j12|cm_d61_j13"
PIN_INPUTS="$(git -C "$R" show "HEAD:$S/$DIR/INPUT_CHECKSUMS.sha256" 2>/dev/null | sha256sum | cut -c1-64)"
PIN_PREREG=0fca97b093f9b987e7a13c589df946f2c134f1c7d466fc7c3e29263508e780c7
PRE=STAGE6_DIAGNOSTIC61_20261125.txt; RESPAT=STAGE6_DIAGNOSTIC61_RESULT; TSV=d61.tsv; EXTRA="d61_functions.tsv"
READER=cs6_d61.py; NEEDS="cs6_d61.py cs6_d61_prepare.py cs6_d61_harvest.sh cs6_d59_prepare.py cs6_d53.py cs6_d56.py cs6_d40.py cs6_d33b.py cs6_audit.py cs6_d47_prepare.py"
TITLE="[retained-final] d61 result: D59's test of the functions at the positive charges of the ten survivors, made again with each ghost 0.01 A beyond its charge"
ANYWAY="${2:-}"
ARC="hpc_all_files_$DIR.tar.xz"
PREREG="$R/$S/$PRE"
WORK="${CS6_WORK:-$HOME/Downloads/cs_${DIR}_harvest_work}"
STATE="$WORK/state.env"
MSG="$WORK/commit_message_results.txt"

say()  { echo "   $*"; }
step() { echo; echo "== $*"; }
die()  { echo; echo "STOPPED: $*"; echo "Nothing has been committed and nothing has been pushed. Nothing was deleted, here or on hpc1."; exit 1; }
sha()  { sha256sum "$1" | cut -c1-64; }
clean() { [ -z "$(git -C "$1" status --porcelain --untracked-files=all)" ]; }
setaside() { local t n=0; t="$2.$(date +%Y%m%d_%H%M%S)"; while [ -e "$t" ]; do n=$((n + 1)); t="$2.$(date +%Y%m%d_%H%M%S).$n"; done; mv "$1" "$t"; }
between() { awk -v a="$1" -v b="$2" 'index($0,a)==1{p=1;next} p&&index($0,b)==1{exit} p&&!/^=+$/&&NF{print}' "$3" | sed 's/^   //'; }

collect() {
  step "1. This PC ($DIR)"
  for t in git ssh scp tar xz sha256sum python3 awk df; do command -v "$t" > /dev/null || die "$t is not installed on this PC."; done
  [ -d "$R/.git" ] || die "$R is not a git repository."
  [ -e "$C/.git" ] || die "$C is not a git repository."
  [ -d "$A" ] || die "the archive folder $A does not exist. Tell Claude where the PC archive lives."
  clean "$R" || die "the results repository has uncommitted or untracked files. If a previous harvest stopped half way, send Claude its output."
  clean "$C" || die "the code repository has uncommitted or untracked files."
  git -C "$R" fetch -q origin || die "could not fetch the results repository."
  git -C "$C" fetch -q origin || die "could not fetch the code repository."
  [ "$(git -C "$R" rev-parse HEAD)" = "$(git -C "$R" rev-parse origin/main)" ] || die "the results repository is not level with GitHub. Pull or push first."
  [ "$(git -C "$C" rev-parse HEAD)" = "$(git -C "$C" rev-parse origin/tier1-realism)" ] || die "the code repository is not level with GitHub. Pull or push first."
  [ -f "$R/$S/$DIR/INPUT_CHECKSUMS.sha256" ] && [ "$(sha "$R/$S/$DIR/INPUT_CHECKSUMS.sha256")" = "$PIN_INPUTS" ] || die "the committed $DIR inputs are not the ones this harvest was written for."
  [ -f "$PREREG" ] && [ "$(sha "$PREREG")" = "$PIN_PREREG" ] || die "the committed $PRE is not the one this harvest was written for."
  for f in $NEEDS; do
    git -C "$C" cat-file -e "HEAD:phase2.2/$f" 2> /dev/null || die "$f is not committed in the code repository. The reader must be committed before it is run on any output."
    git -C "$C" diff --quiet HEAD -- "phase2.2/$f" || die "$f differs from its committed version."
  done
  [ -z "$(find "$R/$S/$DIR" -name job.out -print -quit)" ] || die "the results repository already holds $DIR outputs. This harvest has been run before."
  [ -z "$(find "$R/$S" -maxdepth 1 -name "${RESPAT}_*.txt" -print -quit)" ] || die "a result document for $DIR already exists."
  [ -n "$(find "$R/$S" -maxdepth 1 -name "STAGE6_DIAGNOSTIC59_RESULT_*.txt" -print -quit)" ] || die "D59's result is not committed yet. D61 reads D59's runs N and GN: harvest D59 first."
  [ ! -e "$A/$S/$DIR/$ARC" ] || die "$A/$S/$DIR/$ARC already exists. This harvest has been run before."
  HR="$(git -C "$R" rev-parse HEAD)"; HC="$(git -C "$C" rev-parse HEAD)"
  say "results $HR"; say "code    $HC"; say "archive folder $A"
  if [ -e "$WORK" ]; then setaside "$WORK" "$WORK.old" || die "could not set the old work folder $WORK aside."; fi
  mkdir -p "$WORK/x" || die "could not create $WORK"

  step "2. Have the jobs ended? (plain qstat output from hpc1 follows)"
  Q="$(ssh "$H" "bash -lc 'qstat -u $U'" 2>&1)" || die "the queue on hpc1 could not be read. A failed query is never taken to mean a job has ended."
  echo "$Q" | sed 's/^/   | /'
  [ -n "$Q" ] || say "| (qstat printed nothing: no jobs of yours are queued or running)"
  echo "$Q" | grep -Eqw "$JOBS" && die "a job of D61 is still in the queue. That is expected while they run: run this again when they have all gone."
  ssh "$H" "[ -d $W/$DIR ]" || die "on hpc1, $DIR is missing. The inputs were never uploaded; send Claude the qstat output above."
  SHORT=""
  for L in $LEVELS; do
    DONE="$(ssh "$H" "if [ -s $W/$DIR/run_$L.pbs.out ]; then grep -c '^finished: ' $W/$DIR/run_$L.pbs.out; else echo NOLOG; fi; [ -e $W/$DIR/.running_$L ] && echo LOCK || echo NOLOCK" 2>&1 | tr '\n' ' ')"
    case "$DONE" in
      "1 NOLOCK ") say "$L: the log reaches its last line and the lock is gone: the job ran to its end" ;;
      "NOLOG NOLOCK ") say "$L: there is no log, so the job was never submitted, or has not written its log yet."; SHORT="$SHORT $L" ;;
      *)
        case "$DONE" in *" LOCK "*) say "$L: the .running_$L lock is still there, though the job is not in the queue." ;; esac
        case "$DONE" in "0 "*) say "$L: the log does not reach its last line, which the job script writes when it ends." ;; esac
        case "$DONE" in "NOLOG "*) say "$L: there is no log." ;; esac
        say "$L: so the job was cut short, most likely by its walltime."
        SHORT="$SHORT $L" ;;
    esac
  done
  if [ -n "$SHORT" ]; then
    if [ "$ANYWAY" != "anyway" ]; then
      die "these jobs did not run to their end:$SHORT. If a job was never submitted, run block P5 again: it submits only the jobs that are not yet submitted. If a job was cut short, send Claude this printout: no job script makes a run a second time, so the jobs can be collected as they stand, with block H1A."
    fi
    say "collect anyway: harvesting as it stands. Runs that never ran will be NOT MEASURED."
  fi

  step "3. On hpc1: the inputs that ran, the manifest and the archive (this can take several minutes)"
  cat > "$WORK/remote.sh" <<EOF
set -euo pipefail
cd $W/$DIR
sha256sum -c --quiet INPUT_CHECKSUMS.sha256
echo "INPUTS \$(sha256sum INPUT_CHECKSUMS.sha256 | cut -c1-64)"
if [ ! -e HPC_MANIFEST.sha256 ]; then
  find . -type f ! -name HPC_MANIFEST.sha256 ! -name '*.part' | LC_ALL=C sort | sed 's#^\./##' | xargs sha256sum > HPC_MANIFEST.sha256.part
  mv HPC_MANIFEST.sha256.part HPC_MANIFEST.sha256
fi
sha256sum -c --quiet HPC_MANIFEST.sha256
[ "\$(find . -type f ! -name HPC_MANIFEST.sha256 ! -name '*.part' | wc -l)" = "\$(wc -l < HPC_MANIFEST.sha256)" ]
echo "FILES \$(wc -l < HPC_MANIFEST.sha256)"
echo "BYTES \$(du -sb . | cut -f1)"
cd $W
if [ ! -e $ARC ]; then
  tar -cf - $DIR | xz -T1 -c > $ARC.part
  mv $ARC.part $ARC
fi
echo "ARCHIVE \$(sha256sum $ARC | cut -c1-64)"
echo "MEMBERS \$(tar -tJf $ARC | grep -vc '/\$')"
EOF
  OUT="$(ssh "$H" bash -s < "$WORK/remote.sh" 2>&1)" || { echo "$OUT" | sed 's/^/   | /'; die "the collection on hpc1 failed. If it names a file that did not match, an input or an output changed on hpc1 after it was written."; }
  RIN="$(echo "$OUT" | awk '$1=="INPUTS"{print $2}')"; NFILES="$(echo "$OUT" | awk '$1=="FILES"{print $2}')"
  RSHA="$(echo "$OUT" | awk '$1=="ARCHIVE"{print $2}')"; NMEM="$(echo "$OUT" | awk '$1=="MEMBERS"{print $2}')"; NBYTES="$(echo "$OUT" | awk '$1=="BYTES"{print $2}')"
  [ "$RIN" = "$PIN_INPUTS" ] || die "hpc1 ran inputs whose checksum list is not the committed one."
  [ -n "$NFILES" ] && [ -n "$RSHA" ] && [ -n "$NMEM" ] && [ -n "$NBYTES" ] || die "hpc1's answer could not be read: $OUT"
  [ "$NMEM" = "$((NFILES + 1))" ] || die "the archive on hpc1 holds $NMEM files but the manifest lists $NFILES plus itself."
  say "inputs on hpc1 are the committed ones; $NFILES files in the manifest; archive of $NMEM members"
  NEED=$(( NBYTES / 1000000 + 500 ))
  mkdir -p "$A/$S/$DIR" || die "could not create $A/$S/$DIR"
  FA="$(df --output=avail -BM "$A/$S/$DIR" | tail -n 1 | tr -dc 0-9)"; FW="$(df --output=avail -BM "$WORK" | tail -n 1 | tr -dc 0-9)"
  say "$DIR takes $(( NBYTES / 1000000 )) MB on hpc1. Room needed: $NEED MB in the archive folder (it has $FA MB) and $NEED MB in the work folder $WORK (it has $FW MB)"
  [ -n "$FA" ] && [ "$FA" -ge "$NEED" ] || die "the archive folder has not enough room for the archive. Nothing was downloaded."
  [ -n "$FW" ] && [ "$FW" -ge "$NEED" ] || die "the work folder has not enough room for the unpacked check copy. Nothing was downloaded. Tell Claude: the work folder can be put on the hard drive."

  step "4. Download, and check the copy against hpc1's"
  mkdir -p "$A/$S/$DIR" || die "could not create $A/$S/$DIR"
  [ ! -e "$A/$S/$DIR/$ARC.part" ] || setaside "$A/$S/$DIR/$ARC.part" "$A/$S/$DIR/$ARC.part.old" || die "could not set an old partial download aside."
  scp -q "$H:$W/$ARC" "$A/$S/$DIR/$ARC.part" || die "the download failed. Nothing was placed."
  if [ "$(sha "$A/$S/$DIR/$ARC.part")" != "$RSHA" ]; then
    setaside "$A/$S/$DIR/$ARC.part" "$A/$S/$DIR/$ARC.bad" || die "the downloaded archive does not match hpc1's by sha256, and the bad copy could not be set aside."
    die "the downloaded archive does not match hpc1's by sha256. The bad copy was set aside in $A/$S/$DIR; run this again."
  fi
  mv "$A/$S/$DIR/$ARC.part" "$A/$S/$DIR/$ARC"
  ABYTES="$(wc -c < "$A/$S/$DIR/$ARC" | tr -d ' ')"
  say "$A/$S/$DIR/$ARC  $ABYTES bytes  sha256 $RSHA"

  step "5. Every member of the archive against hpc1's manifest, and against the repository's committed files"
  tar -xJf "$A/$S/$DIR/$ARC" -C "$WORK/x" || die "the archive could not be unpacked."
  ( cd "$WORK/x/$DIR" && sha256sum -c --quiet HPC_MANIFEST.sha256 ) || die "a member of the archive does not match hpc1's manifest."
  [ "$(find "$WORK/x/$DIR" -type f | wc -l | tr -d ' ')" = "$NMEM" ] || die "the unpacked archive does not hold $NMEM files."
  NTRK=0
  while IFS= read -r f; do
    cmp -s "$WORK/x/$DIR/${f#"$S/$DIR/"}" "$R/$f" || die "${f#"$S/$DIR/"} on hpc1 is not the committed file."
    NTRK=$((NTRK + 1))
  done < <(git -C "$R" ls-files "$S/$DIR")
  [ "$NTRK" -gt "$NIN" ] || die "the repository tracks $NTRK files under $DIR; more than $NIN were expected."
  say "all $NMEM members match the manifest; all $NTRK committed files of $DIR (the $NIN checksummed inputs among them) are the ones hpc1 ran"

  step "6. Place the text outputs and write the records"
  python3 - "$WORK/x/$DIR" "$R" "$S/$DIR" "$A" "$S/$DIR/$ARC" "$ABYTES" "$RSHA" "$NMEM" "$HR" "$ARC" "$NRUN" <<'PY' || die "placing the outputs failed; see the line above. Send Claude this output before running anything else."
import hashlib, os, shutil, subprocess, sys
src, repo, rel, arch, arcrel, abytes, asha, nmem, head, arcname, nrun = sys.argv[1:12]
dst = os.path.join(repo, rel)
name = os.path.basename(rel)
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 16), b""):
            h.update(b)
    return h.hexdigest()
tracked = subprocess.run(["git", "-C", repo, "ls-files", rel], stdout=subprocess.PIPE, check=True).stdout.decode().split("\n")
committed = set(os.path.relpath(t, rel) for t in tracked if t)
files = sorted(os.path.relpath(os.path.join(d, f), src) for d, _, fs in os.walk(src) for f in fs)
LEVELS = ("J01", "J02", "J03", "J04", "J05", "J06", "J07", "J08", "J09", "J10", "J11", "J12", "J13")
TOP = tuple("run_%s.pbs.out" % l for l in LEVELS) + ("HPC_MANIFEST.sha256",)
def earlier_log(p):
    return any(p.startswith("run_%s.pbs.out.part" % l) and p[len("run_%s.pbs.out.part" % l):].isdigit() for l in LEVELS)
PER = ("job.out", "job.property.txt")
inv, new = [], []
for p in files:
    parts = p.split("/")
    if p in committed:
        cat = "in_repo"
    elif (len(parts) == 2 and parts[1] in PER) or p in TOP or earlier_log(p):
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
if len(runs) != int(nrun):
    sys.exit("   REFUSED: the repository holds %d run directories under %s, not %s." % (len(runs), name, nrun))
outs = [r for r in runs if os.path.isfile(os.path.join(dst, r, "job.out"))]
norm = [r for r in outs if "ORCA TERMINATED NORMALLY" in open(os.path.join(dst, r, "job.out"), errors="replace").read()]
pc = sorted(set(os.path.basename(p) for b, c, p in inv if c == "pc_only"))
with open(os.path.join(dst, "README.txt"), "w") as f:
    f.write("""%s: the %d run(s) kept as physical files, then read (cs6_d61_harvest.sh).

Rule (Jessica, 18 October 2026): every output file of a finished job is downloaded and kept as a
physical file, not only reported.

1. PC archive (PC_ARCHIVES.tsv), in the archive folder on Jessica's PC: all %d files of %s as they
   stood on hpc1 after the jobs had ended, together with hpc1's own HPC_MANIFEST.sha256, which is
   %s members in one %s. hpc1 wrote HPC_MANIFEST.sha256 over its own files
   before the download; the archive was checked against hpc1's copy by sha256 after the download,
   and every archive member was checked against the manifest again where the archive lies.
2. In this repo: job.out and job.property.txt for each run that wrote them; the 13 PBS
   logs
   run_<job>.pbs.out, and run_<job>.pbs.out.part<N> for each earlier submission if a job
   was submitted more than once; HPC_MANIFEST.sha256, INPUT_CHECKSUMS.sha256,
   INVENTORY.tsv (every file with its size and category) and OUTPUTS.sha256.
3. NOT in this repo: %s.
   Wavefunction and scratch files are excluded by this repository's own .gitignore, the split
   decided on 28 October 2026 for Diagnostic 14: these are fixed-geometry single points whose
   reading uses the text files kept here, and the inputs that would regenerate them are committed beside
   them. Every one of these files is in the PC archive with its sha256 in HPC_MANIFEST.sha256.
4. The harvest deletes nothing on hpc1.

job.out is present for %d of the %d run(s) and reports normal termination in %d. A run that did
not finish is a result and is kept as it is; the reader reports it.

The inputs were committed and pushed before any job was submitted; the harvest started from
results %s. Every committed file hpc1 ran was checked against the committed one twice: by
sha256sum -c on hpc1 before the archive was made, and byte for byte between the archive's copies
and this repository's.

The archive is named %s.
""" % (name, len(runs), len(files) - 1, name, nmem, arcname, ", ".join(pc) if pc else "nothing", len(outs), len(runs), len(norm), head[:8], arcname))
print("   placed %d files; job.out present for %d of %d run(s), normal termination in %d" % (len(new), len(outs), len(runs), len(norm)))
print("   kept in the PC archive only: %s" % (", ".join(pc) if pc else "nothing"))
PY

  step "7. Read: the committed reader, $READER"
  RES="${RESPAT}_$(date +%Y%m%d).txt"
  if ! python3 "$C/phase2.2/$READER" --res "$R" --code "$C" --prereg "$PREREG" --result "$R/$S/$RES" --tsv "$R/$S/$TSV" > "$WORK/reader.out" 2>&1; then
    cat "$WORK/reader.out"
    echo
    echo "THE READER STOPPED AT A GATE, so there is no result document."
    echo "The archive is on your PC and the outputs are in place in the results repository, uncommitted."
    echo "Do NOT run the commit step. Send Claude everything printed above."
    exit 2
  fi
  cat "$WORK/reader.out"
  [ -s "$R/$S/$RES" ] && [ -s "$R/$S/$TSV" ] && [ -s "$R/$S/d61_functions.tsv" ] || die "the reader ended without writing its result."

  {
    echo "$TITLE"
    echo
    echo "Read by the committed cs6_d61.py from $RES, which holds every number below. D59's R59.1 made again (D59 NOTE 1): the functions of D50 at the positive charges of the ten survivors of D56 with the largest lowering, each ghost centre 0.01 A beyond its charge, read as D53 read its arrangements at CAM-B3LYP-D3BJ and wB97X-V (def2-TZVPD) with Hartree-Fock beside, F and GF from d61/ and N and GN from D59's committed runs; in vacuo."
    echo
    grep -E "WITH FUNCTIONS AT ITS POSITIVE CHARGES|NOT SETTLED WITH FUNCTIONS" "$R/$S/$RES" | sed 's/^ *//'
    echo
    echo "Reported beside D56's word and D59's walled word, and nothing more: no design is named, kept or dropped on these readings and no word of D56 or D59 changes. Whether a design that takes electrons may be quoted is for the amendment that follows this result (A11.5 item 5)."
    echo
    awk -v RS= '{p=$0} END{gsub(/\n/," ",p); print p}' "$R/$S/$RES"
    echo
    echo "Kept as physical files: the text outputs of the runs here, and all $NMEM files in the PC archive (PC_ARCHIVES.tsv, sha256 $RSHA). Nothing was deleted on hpc1."
    echo
    echo "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
    echo "Claude-Session: https://claude.ai/code/session_01F3vG9e2HqrawvyQuL6J8Ke"
  } > "$MSG"
  {
    echo "HR=$HR"; echo "HC=$HC"; echo "RES=$RES"; echo "RSHA=$RSHA"; echo "NMEM=$NMEM"
  } > "$STATE"

  step "COLLECTED ($DIR). Nothing is committed yet."
  say "The result is $R/$S/$RES"
  say "The commit message that the next step will use is $MSG"
  say "The next step is: bash $C/phase2.2/cs6_d61_harvest.sh commit"
}

commit() {
  [ -s "$STATE" ] && [ -s "$MSG" ] || die "collect has not finished successfully for $DIR, so there is nothing to commit."
  # shellcheck disable=SC1090
  . "$STATE"
  [ "$(git -C "$R" rev-parse HEAD)" = "$HR" ] || die "the results repository has moved since collect ran."
  [ "$(git -C "$C" rev-parse HEAD)" = "$HC" ] || die "the code repository has moved since collect ran."
  clean "$C" || die "the code repository has uncommitted or untracked files."
  [ -s "$R/$S/$RES" ] || die "the result document is missing."
  BAD="$(git -C "$R" status --porcelain --untracked-files=all | grep -vc -E "^\?\? $S/($DIR/|$TSV\$|$RES\$|d61_functions\.tsv\$)" || true)"
  [ "$BAD" = 0 ] || die "the results repository holds changes that are not part of this harvest."
  [ -n "$(git -C "$R" status --porcelain --untracked-files=all)" ] || die "there is nothing to commit in the results repository."
  ( cd "$R/$S/$DIR" && sha256sum -c --quiet OUTPUTS.sha256 ) || die "a placed output no longer matches OUTPUTS.sha256."
  ADD=""; for x in $EXTRA; do [ -e "$R/$S/$x" ] && ADD="$ADD $S/$x"; done
  git -C "$R" add "$S/$DIR" "$S/$TSV" "$S/$RES" $ADD || die "git add failed in the results repository."
  git -C "$R" commit -q -F "$MSG" || die "the results commit failed."
  clean "$R" || die "something is left uncommitted after the commit. Do not push; send Claude the output of git status."
  echo
  git -C "$R" log -1 --format='results committed: %h %s'
  echo "The code repository is unchanged: the reader and this script were committed before any output was read."
  echo "Nothing has been pushed."
}

case "${1:-}" in
  collect) collect ;;
  commit)  commit ;;
  *) echo "usage: bash cs6_d61_harvest.sh <collect|collect anyway|commit>"; exit 1 ;;
esac
