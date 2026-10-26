#!/usr/bin/env python3
"""Stage 6 DIAGNOSTIC 8: the wavefunction-health audit re-run over the corpus it should have had.

NO NEW CALCULATION AND NO NETWORK. Diagnostic 6 is not edited and DIAGNOSTIC6_RESULT.txt is not
touched. This is a separate audit, with a separate result file, that applies D6.3's classification
rule UNCHANGED to a corrected corpus. STAGE6_ERRATUM4_20261026.txt names four defects in
cs6_audit.py; this script fixes them in three nested stages so that the direction of each fix is
visible separately, because two of them pull opposite ways.

  stage 1   D6 exactly as published. Not a new reading: a GATE.
  stage 2   + controls restricted to charge-free SINGLE POINTS, and the four run categories made
            disjoint.                                                             (E4.4 i and ii)
  stage 3   + loose .xz outputs walked, and point charges detected from Ne> centres as well as
            from a point-charge file.                                             (E4.5 and E4.6)

Stage 3's two fixes MUST be made together, which is Erratum 4 E4.6's binding constraint: the
walker alone adds runs whose charges are Ne> centres, and without the charge fix those runs are
read as charge-free CONTROLS, which would make the audit more permissive rather than more complete.

THE GATE. Before any new reading is reported, this script reconstructs the ENTIRE text of
DIAGNOSTIC6_RESULT.txt from stage 1 and compares it byte for byte with the committed file. The
corpus reader, the control medians and the classification are therefore all validated against the
published artefact itself, not against any number retyped from it. If the reconstruction differs,
THAT IS THE FINDING: the script stops and reports the first differing line. Nothing else is
printed, because nothing else would be interpretable.

cs6_audit.read() and every one of D6's thresholds are IMPORTED from cs6_audit.py and used
unchanged. Nothing is retyped. verdict() cannot be imported - it is a closure inside cs6_audit's
main() - so it is reproduced here line for line against D6's imported constants, and the gate
above is what proves the reproduction is faithful.

usage: cs8_audit.py <results-repo> <code-repo> <out-dir>
"""
import sys, re, tarfile, lzma, math, hashlib, statistics as st
from pathlib import Path

# importing cs6_audit would otherwise drop a __pycache__ directory into the code repo and leave it
# dirty, which is both a side effect on a repo this script is supposed to only read and a thing that
# would make the next run refuse. Set before the import, not after.
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cs6_audit as D6                      # read(), D_LW, R_DIP, GAP_ABS, R_GAP, R_CYC, MM_NPC

KE = 14.39964548          # e^2 / (4 pi eps0), in eV * Angstrom; E = KE*q/r^2 gives V/Angstrom
ANG = "CARTESIAN COORDINATES (ANGSTROEM)"

# Erratum 4 E4.6 quotes the Ne> centres of L0_E5_R verbatim; the three positions matched
# L1_E5_R/combined.xyz to eight decimals, which is where this pattern comes from.
RE_NE = re.compile(r"^\s*Ne>\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s*$")
RE_AT = re.compile(r"^\s*([A-Za-z][A-Za-z]?)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s*$")


def stop(m):
    print("STOP: %s" % m, file=sys.stderr)
    raise SystemExit(1)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# --------------------------------------------------------------------- the two added readings
def ne_centres(s):
    """Every Ne> point charge in the LAST Angstroem coordinate block, and in the file as a whole.

    Returns (block_atoms, block_charges, n_file). A charge is counted only if it is non-zero: a
    zero-charge Ne> centre is a ghost centre and is not a designed charge.
    """
    seen, nfile = set(), 0
    for l in s.splitlines():
        m = RE_NE.match(l)
        if m:
            q, x, y, z = (float(g) for g in m.groups())
            if abs(q) > 0.0:
                k = (round(q, 6), round(x, 6), round(y, 6), round(z, 6))
                if k not in seen:
                    seen.add(k)
                    nfile += 1
    k = s.rfind(ANG)
    atoms, chg = [], []
    if k >= 0:
        for l in s[k:].splitlines()[2:]:
            if not l.strip():
                if atoms or chg:
                    break
                continue
            m = RE_NE.match(l)
            if m:
                q, x, y, z = (float(g) for g in m.groups())
                if abs(q) > 0.0:
                    chg.append((q, x, y, z))
                continue
            m = RE_AT.match(l)
            if m:
                atoms.append((float(m.group(2)), float(m.group(3)), float(m.group(4))))
                continue
            break
    return atoms, chg, nfile


def field(atoms, chg):
    """Peak |E| over the real atoms, in V/Angstrom, and the closest charge-to-atom approach."""
    if not atoms or not chg:
        return None, None
    best, near = 0.0, 1e9
    for ax, ay, az in atoms:
        ex = ey = ez = 0.0
        for q, cx, cy, cz in chg:
            dx, dy, dz = ax - cx, ay - cy, az - cz
            r2 = dx * dx + dy * dy + dz * dz
            if r2 < 1e-10:
                continue
            r = math.sqrt(r2)
            near = min(near, r)
            f = KE * q / (r2 * r)
            ex += f * dx; ey += f * dy; ez += f * dz
        best = max(best, math.sqrt(ex * ex + ey * ey + ez * ez))
    return best, (None if near > 1e8 else near)


# --------------------------------------------------------------------- the corpus
def pc_beside(pcmap, outs, name):
    """The .pc file belonging to this output inside the same archive.

    A charge vector attributed to the wrong run is worse than none at all, so the match must be
    unambiguous: the same stem, or - the ORCA convention of one job per directory - a directory
    holding exactly one output and exactly one .pc file. Anything else returns nothing.
    """
    stem = name.rsplit(".", 1)[0]
    if stem + ".pc" in pcmap:
        return pcmap[stem + ".pc"]
    dr = lambda k: k.rsplit("/", 1)[0] if "/" in k else ""
    d = dr(name)
    same = [v for k, v in pcmap.items() if dr(k) == d]
    if len(same) == 1 and sum(1 for k in outs if dr(k) == d) == 1:
        return same[0]
    return None


def pc_on_disk(path):
    """The .pc file sitting next to a loose output, plain or xz."""
    stem = path.name
    for suf in (".out.xz", ".xz"):
        if stem.endswith(suf):
            stem = stem[:-len(suf)]
            break
    for cand in (path.parent / (stem + ".pc"), path.parent / (stem + ".pc.xz")):
        if cand.is_file():
            try:
                if cand.suffix == ".xz":
                    with lzma.open(str(cand)) as fh:
                        return fh.read().decode(errors="replace")
                return cand.read_text(errors="replace")
            except Exception:
                return None
    return None


def read_pc(text):
    """An ORCA point-charge file: a count on the first line, then q x y z per line."""
    L = [l.split() for l in text.splitlines() if l.strip()]
    if not L or len(L[0]) != 1:
        return []
    try:
        n = int(L[0][0])
    except ValueError:
        return []
    out = []
    for f in L[1:]:
        if len(f) < 4:
            continue
        try:
            out.append((float(f[0]), float(f[1]), float(f[2]), float(f[3])))
        except ValueError:
            return []
    return out if len(out) == n else []


def collect(root):
    """One pass. Every record is tagged with how it was found, so a stage can filter on it."""
    root = Path(root)
    if not (root / ".git").is_dir():
        stop("%s is not a git repo" % root)
    R, loose_out = [], 0

    def take(text, arc, run, src, pcf=None):
        d = D6.read(text)
        if not d:
            return
        at, ch, nf = ne_centres(text)
        d["arc"], d["run"], d["src"] = arc, run, src
        d["nne_blk"], d["nne"] = len(ch), max(nf, len(ch))
        d["chg"] = ch if len(ch) == d["nne"] else []
        d["qsrc"] = "Ne>" if d["chg"] else "-"
        # a point-charge FILE committed beside the output is the charges themselves, not an echo
        if not d["chg"] and d["npc"] and pcf:
            pc = read_pc(pcf)
            if len(pc) == d["npc"]:
                d["chg"], d["qsrc"] = pc, ".pc"
        d["geo"] = hashlib.sha256(
            ("|".join("%.4f,%.4f,%.4f" % t for t in at)).encode()).hexdigest()[:12] if at else "-"
        if d["chg"]:
            d["qid"] = hashlib.sha256(
                ("|".join("%.4f,%.4f,%.4f,%.4f" % t for t in sorted(d["chg"]))).encode()).hexdigest()[:12]
            d["sumq"] = sum(abs(q) for q, _, _, _ in d["chg"])
            d["maxq"] = max(abs(q) for q, _, _, _ in d["chg"])
            d["sq2"] = sum(q * q for q, _, _, _ in d["chg"])
            d["Epk"], d["rmin"] = field(at, d["chg"])
            d["nat"] = len(at)
        else:
            d["qid"] = None
            d["sumq"] = d["maxq"] = d["sq2"] = d["Epk"] = d["rmin"] = None
            d["nat"] = len(at)
        R.append(d)

    for a in sorted(p for p in root.rglob("*.tar.xz") if ".git" not in p.parts):
        try:
            t = tarfile.open(str(a))
        except Exception:
            continue
        pcmap, outs = {}, []
        for m in t.getmembers():
            if not m.isfile():
                continue
            if m.name.endswith(".pc"):
                try:
                    pcmap[m.name] = t.extractfile(m).read().decode(errors="replace")
                except Exception:
                    pass
            elif m.name.endswith(".out") and not m.name.endswith(".pbs.out"):
                outs.append(m.name)
        for m in t.getmembers():
            if not m.isfile() or not m.name.endswith(".out") or m.name.endswith(".pbs.out"):
                continue
            try:
                s = t.extractfile(m).read().decode(errors="replace")
            except Exception:
                continue
            take(s, str(a.relative_to(root)), m.name, "archive", pc_beside(pcmap, outs, m.name))

    # E4.5: completed runs committed as loose .xz were never opened at all.
    for a in sorted(p for p in root.rglob("*.xz")
                    if ".git" not in p.parts and not p.name.endswith(".tar.xz")):
        try:
            with lzma.open(str(a)) as fh:
                s = fh.read().decode(errors="replace")
        except Exception:
            continue
        if s.lstrip().startswith("\x1f") or "FINAL SINGLE POINT ENERGY" not in s:
            continue
        n0 = len(R)
        take(s, str(a.relative_to(root)), a.name[:-3], "loose", pc_on_disk(a))
        loose_out += len(R) - n0

    # counted and reported, never classified: whether to include them is a separate decision
    plain = 0
    for p in root.rglob("*.out"):
        if ".git" in p.parts or p.name.endswith(".pbs.out"):
            continue
        try:
            if D6.read(p.read_text(errors="replace")):
                plain += 1
        except Exception:
            continue
    return R, loose_out, plain


# --------------------------------------------------------------------- D6.3, reproduced exactly
def controls(good, strict):
    CTL = {}
    for nb in {d["nb"] for d in good}:
        g = [d for d in good if d["nb"] == nb and d["pc"] == 0 and (d["nopt"] == 0 if strict else True)]
        if g:
            CTL[nb] = {k: st.median([x[k] for x in g if x[k] is not None])
                       for k in ("cyc", "gap", "maxlw", "dip") if any(x[k] is not None for x in g)}
            CTL[nb]["n"] = len(g)
    return CTL


def verdict(d, CTL):
    c = CTL.get(d["nb"])
    if not c or "gap" not in c:
        return "NO CONTROL", []
    moved, sick = [], []
    if d["maxlw"] is not None and "maxlw" in c and d["maxlw"] - c["maxlw"] >= D6.D_LW:
        moved.append("Loewdin +%.3f e over control" % (d["maxlw"] - c["maxlw"]))
    if d["dip"] is not None and c.get("dip") and d["dip"] / c["dip"] >= D6.R_DIP:
        moved.append("dipole %.1fx control" % (d["dip"] / c["dip"]))
    if d["gap"] is not None and d["gap"] < D6.GAP_ABS and d["gap"] < D6.R_GAP * c["gap"]:
        sick.append("gap %.3f eV" % d["gap"])
    if d["cyc"] is not None and c.get("cyc") and d["cyc"] >= D6.R_CYC * c["cyc"]:
        sick.append("%d SCF cycles, %.1fx control" % (d["cyc"], d["cyc"] / c["cyc"]))
    if moved and sick:
        return "COLLAPSED", moved + sick
    if sick:
        return "STRAINED", sick
    if moved:
        return "POLARISED", moved
    return "CLEAN", []


def stage(R, n):
    """Apply stage n's corpus, charge rule, control filter and categories. Returns a dict."""
    # each stage works on its own copies, so one stage can never overwrite another's verdicts
    pool = [dict(d) for d in R if n >= 3 or d["src"] == "archive"]
    for d in pool:
        d["pc"] = max(d["npc"], d["nne"]) if n >= 3 else d["npc"]
    good = [d for d in pool if d["norm"] and not d["notconv"]]
    CTL = controls(good, strict=(n >= 2))
    for d in good:
        d["v"], d["w"] = verdict(d, CTL)
    if n == 1:                                           # D6's own categories, overlaps and all
        cf = [d for d in good if d["pc"] == 0]
        sp = [d for d in good if 0 < d["pc"] <= D6.MM_NPC and d["nopt"] == 0]
        mm = [d for d in good if d["pc"] > D6.MM_NPC]
        go = [d for d in good if d["pc"] > 0 and d["nopt"] > 0]
        cfo = []
    else:                                                # E4.4(ii): five disjoint classes
        cf = [d for d in good if d["pc"] == 0 and d["nopt"] == 0]
        cfo = [d for d in good if d["pc"] == 0 and d["nopt"] > 0]
        sp = [d for d in good if 0 < d["pc"] <= D6.MM_NPC and d["nopt"] == 0]
        go = [d for d in good if 0 < d["pc"] <= D6.MM_NPC and d["nopt"] > 0]
        mm = [d for d in good if d["pc"] > D6.MM_NPC]
    return dict(n=n, pool=pool, good=good, CTL=CTL, cf=cf, cfo=cfo, sp=sp, mm=mm, go=go)


# --------------------------------------------------------------------- the gate
def render_d6(S):
    """DIAGNOSTIC6_RESULT.txt, rebuilt from stage 1. Copied from cs6_audit.main() so that the
    comparison is total: every number D6 printed has to come back out of this code."""
    R, good, CTL = S["pool"], S["good"], S["CTL"]
    cf, sp, mm, go = S["cf"], S["sp"], S["mm"], S["go"]
    O = ["STAGE 6 - DIAGNOSTIC 6 RESULT: a wavefunction-health audit of the whole committed charge-design record",
         "",
         "Measured by cs6_audit.py. NO NEW CALCULATION. Four fields are read out of every committed ORCA output",
         "on this dianion that Diagnostics 1 to 5 never read: the SCF cycle count, the HOMO-LUMO gap, the dipole",
         "magnitude and the largest Loewdin charge on a light atom. The rule is STAGE6_DIAGNOSTIC6_20261026.txt",
         "D6.3, applied unchanged.", "",
         "   %d completed single points and optimisations on this dianion (118 electrons, 59 occupied)" % len(good),
         "   %d carry no point charges and are the controls; %d are point-charge single points; %d are"
         % (len(cf), len(sp), len(go)),
         "   point-charge geometry optimisations, reported separately because a charge-free control is NOT",
         "   like-for-like once the geometry has moved. A further %d runs carry more than %d point charges and" % (len(mm), D6.MM_NPC),
         "   are QM/MM calculations whose charges are an MM environment rather than a designed arrangement;",
         "   they are excluded from every table below.",
         "   %d outputs were skipped as not completed or not this system." % (len(R) - len(good)), "",
         "FOUR WORDS, AND WHAT EACH MEANS", "",
         "   CLEAN      no indicator departs from the same-basis charge-free control.",
         "   POLARISED  the density has moved - Loewdin up at least %.2f e over control, or the dipole at least"
         % D6.D_LW,
         "              %.1f times it - but the SCF is healthy. THIS IS NOT A FAULT. It is what a strong"
         % D6.R_DIP,
         "              external field is supposed to do, and it is reported so that it is not mistaken for one.",
         "   STRAINED   the SCF is unwell - gap below %.1f eV AND below %.0f per cent of its own control's gap,"
         % (D6.GAP_ABS, 100 * D6.R_GAP),
         "              or at least %.0f times the control's cycle count - without demonstrable relocation."
         % D6.R_CYC,
         "              A near-degeneracy: single-reference closed-shell DFT is not reliable there.",
         "   COLLAPSED  both at once. The density has moved AND the SCF is unwell. This is the signature of",
         "              density being captured by a positive point charge that has no Pauli wall for electrons",
         "              (Finding 4), the effect Marefat Khah et al. JCTC 2020 doi:10.1021/acs.jctc.9b01162",
         "              describe as unphysical leaking of the electron density.", "",
         "THE CONTROLS EVERY THRESHOLD IS ANCHORED TO", "",
         "   %5s %6s %5s %10s %9s %9s   %s" % ("nb", "basis", "n", "med cyc", "med gap", "med maxLw", "med dip")]
    for nb in sorted(CTL):
        c = CTL[nb]
        bs = [d["bas"] for d in cf if d["nb"] == nb]
        bs = max(set(bs), key=bs.count) if bs else "?"
        O.append("   %5d %-6s %5d %10.0f %9.3f %9.3f   %6.2f"
                 % (nb, bs, c["n"], c.get("cyc", 0), c.get("gap", 0), c.get("maxlw", 0), c.get("dip", 0)))
    nbad = sum(1 for d in cf if d["v"] != "CLEAN")
    O += ["",
          "   SPECIFICITY CHECK: of the %d charge-free controls, %d are classified anything but CLEAN."
          % (len(cf), nbad),
          "   %s" % ("A control that trips the classifier would mean the thresholds are too tight; none does."
                     if nbad == 0 else "THAT IS NOT ZERO AND THE THRESHOLDS ARE TOO TIGHT - say so in any use of this table."),
          "",
          "POINT-CHARGE SINGLE POINTS, BY BASIS SET", "",
          "   %5s %-12s %6s %8s %7s %9s %8s %9s" % ("nb", "basis", "runs", "npc med", "clean", "polarised", "strained", "collapsed")]
    for nb in sorted({d["nb"] for d in sp}):
        g = [d for d in sp if d["nb"] == nb]
        bs = [d["bas"] for d in g]
        n = lambda v: sum(1 for d in g if d["v"] == v)
        O.append("   %5d %-12s %6d %8d %7d %9d %8d %9d"
                 % (nb, max(set(bs), key=bs.count), len(g), st.median([d["pc"] for d in g]),
                    n("CLEAN"), n("POLARISED"), n("STRAINED"), n("COLLAPSED")))
    O += ["", "AT HOW MANY CHARGES DOES EACH BASIS COLLAPSE?", "",
          "   This is the reading that answers why the representation passed its own tests and then failed.", "",
          "   %5s %-12s %10s %28s" % ("nb", "basis", "collapsed", "number of point charges in them")]
    for nb in sorted({d["nb"] for d in sp if d["v"] == "COLLAPSED"}):
        g = [d for d in sp if d["nb"] == nb and d["v"] == "COLLAPSED"]
        v = sorted(d["pc"] for d in g)
        bs = [d["bas"] for d in g]
        O.append("   %5d %-12s %10d %28s" % (nb, max(set(bs), key=bs.count), len(g),
                 "%d to %d, median %d" % (v[0], v[-1], st.median(v))))
    for nb in sorted({d["nb"] for d in sp}):
        g = [d for d in sp if d["nb"] == nb]
        if not any(d["v"] == "COLLAPSED" for d in g):
            v = sorted(d["pc"] for d in g)
            bs = [d["bas"] for d in g]
            O.append("   %5d %-12s %10d %28s" % (nb, max(set(bs), key=bs.count), 0,
                     "none, over npc %d to %d" % (v[0], v[-1])))
    O += ["", "INDICATORS AGAINST THE NUMBER OF POINT CHARGES, def2-SVP ONLY (nb 264)", "",
          "   def2-SVP has no diffuse functions, so it is the basis every committed Stage 6 design was",
          "   evaluated at, and the one whose quiet failure matters most.", "",
          "   %8s %6s %10s %10s %11s %12s %10s"
          % ("npc", "runs", "med cyc", "med gap", "worst gap", "med maxLw", "collapsed")]
    for lo, hi, lab in ((1, 1, "1"), (2, 4, "2-4"), (5, 9, "5-9"), (10, 13, "10-13"),
                        (14, 29, "14-29"), (30, 10 ** 6, ">=30")):
        g = [d for d in sp if d["nb"] == 264 and lo <= d["pc"] <= hi]
        if not g:
            continue
        md = lambda k: (lambda v: "%.3f" % st.median(v) if v else "-")([x[k] for x in g if x[k] is not None])
        O.append("   %8s %6d %10s %10s %11s %12s %10d"
                 % (lab, len(g), md("cyc"), md("gap"),
                    "%.3f" % min([x["gap"] for x in g if x["gap"] is not None], default=0),
                    md("maxlw"), sum(1 for d in g if d["v"] == "COLLAPSED")))
    O += ["", "EVERY COLLAPSED OR STRAINED SINGLE POINT", "",
          "   %-38s %-28s %4s %4s %4s  %s" % ("archive", "run", "npc", "nb", "cyc", "verdict and why")]
    for d in sorted([x for x in sp if x["v"] in ("COLLAPSED", "STRAINED")],
                    key=lambda x: (x["v"], x["arc"], x["run"])):
        O.append("   %-38s %-28s %4d %4d %4s  %s: %s"
                 % (d["arc"][-38:], d["run"].split("/")[-1][:28], d["pc"], d["nb"],
                    d["cyc"] if d["cyc"] is not None else "-", d["v"], "; ".join(d["w"])))
    O += ["", "THE POINT-CHARGE GEOMETRY OPTIMISATIONS, REPORTED APART", "",
          "   A charge-free control is not like-for-like once the geometry has moved, so these verdicts carry",
          "   much less weight than the single points above and must not be quoted as collapses on their own.", ""]
    for v in ("CLEAN", "POLARISED", "STRAINED", "COLLAPSED", "NO CONTROL"):
        n = sum(1 for d in go if d["v"] == v)
        if n:
            O.append("   %-11s %4d" % (v, n))
    O += ["", "WHAT THIS DOES NOT TEST", "",
          "   Nothing was computed, so this cannot test any hypothesis the committed outputs do not already",
          "   contain. It cannot say what the correct barrier is for a run it calls COLLAPSED or STRAINED - only",
          "   that the number should not be quoted. It does not establish a usable charge model; a Pauli wall",
          "   for electrons is the fix and it is untested here. The four thresholds were chosen AFTER the Stage 6",
          "   outputs had been inspected, which is why they are anchored to in-corpus controls and why the",
          "   control table and the specificity check are printed: the classification is reproducible and",
          "   checkable, but it is not a blind test. A Loewdin charge is basis-dependent and is used here only",
          "   as a difference against a same-basis control, never as an absolute quantity.", ""]
    return "\n".join(O) + "\n"


# --------------------------------------------------------------------- the new readings
def key(d):
    return (d["arc"], d["run"])


def census(O, S, lab):
    O.append("   %-9s %7d %9d %8d %8d %8d %7d"
             % (lab, len(S["good"]), len(S["cf"]), len(S["cfo"]), len(S["sp"]),
                len(S["go"]), len(S["mm"])))


def vcount(g, v):
    return sum(1 for d in g if d["v"] == v)


def render_d8(S1, S2, S3, R, nloose, nplain, d6path, d6sha, cs6sha, capped=160):
    O = ["STAGE 6 - DIAGNOSTIC 8 RESULT: the wavefunction-health audit re-run over the corpus it should have had",
         "",
         "Measured by cs8_audit.py. NO NEW CALCULATION AND NO NETWORK. DIAGNOSTIC6_RESULT.txt is not edited and",
         "STAGE6_DIAGNOSTIC6_20261026.txt is not amended. D6.3's classification rule is applied UNCHANGED: every",
         "threshold is imported from cs6_audit.py, and cs6_audit.read() parses every output here exactly as it",
         "did there. One variable changes between stages, and it is the corpus, never the rule.",
         "",
         "   cs6_audit.py              sha256 %s" % cs6sha,
         "   DIAGNOSTIC6_RESULT.txt    sha256 %s" % d6sha,
         "   read from                 %s" % d6path,
         "",
         "THE GATE, WHICH PASSED",
         "",
         "   Stage 1 below is D6 exactly as published. Before any new reading was computed, the ENTIRE text of",
         "   DIAGNOSTIC6_RESULT.txt was rebuilt from stage 1 and compared with the committed file byte for byte.",
         "   It matched. Every number D6 printed therefore comes back out of this script's corpus reader, its",
         "   control medians and its reproduction of D6.3 - so the stages below differ from D6 only where a",
         "   defect named in Erratum 4 has been corrected, and nowhere else. Had it not matched, this file would",
         "   not exist: the script stops there and reports the first differing line, because no later reading",
         "   would be interpretable.",
         "",
         "WHAT EACH STAGE FIXES",
         "",
         "   stage 1   D6 as published. The gate, not a reading.",
         "   stage 2   + controls restricted to charge-free SINGLE POINTS (D6 pooled charge-free geometry",
         "             optimisations into the control medians, against its own D6.1), and the run categories",
         "             made disjoint (D6 counted a QM/MM optimisation in two categories at once). E4.4 i, ii.",
         "   stage 3   + loose .xz outputs opened (E4.5: cs6_audit.py walked *.tar.xz only), and point charges",
         "             detected from Ne> centres as well as from a point-charge file (E4.6).",
         "",
         "   Stage 3's two fixes are made TOGETHER and cannot be separated. They pull opposite ways: the walker",
         "   ADDS runs, and most of the added runs carry their charges as Ne> centres, so without the charge fix",
         "   they would enter as charge-free CONTROLS and loosen every threshold in the table. Reporting the",
         "   stages nested is what keeps the two directions visible instead of netted.",
         "",
         "(a) THE CORPUS AT EACH STAGE",
         "",
         "   %-9s %7s %9s %8s %8s %8s %7s"
         % ("stage", "kept", "controls", "cf-opt", "pc-sp", "pc-opt", "QM/MM")]
    for S, lab in ((S1, "1 (D6)"), (S2, "2"), (S3, "3")):
        census(O, S, lab)
    tot1 = len(S1["cf"]) + len(S1["cfo"]) + len(S1["sp"]) + len(S1["go"]) + len(S1["mm"])
    tot3 = len(S3["cf"]) + len(S3["cfo"]) + len(S3["sp"]) + len(S3["go"]) + len(S3["mm"])
    O += ["",
          "   Stage 1's five columns sum to %d against %d runs kept: %s. Stage 3's sum to %d against %d: %s."
          % (tot1, len(S1["good"]), "D6's categories overlap" if tot1 != len(S1["good"]) else "no overlap",
             tot3, len(S3["good"]), "STILL OVERLAPPING" if tot3 != len(S3["good"]) else "a partition"),
          "",
          "   controls  charge-free single points (stage 1: charge-free runs of any kind, which is the defect)",
          "   cf-opt    charge-free geometry optimisations, split out at stage 2 and never used as controls",
          "   pc-sp     point-charge single points, 1 to %d charges - the table everything is read off" % D6.MM_NPC,
          "   pc-opt    point-charge geometry optimisations, reported apart",
          "   QM/MM     more than %d charges: an MM environment, not a design" % D6.MM_NPC,
          "",
          "   %d completed runs on this dianion were added by opening loose .xz outputs." % nloose,
          "   %d further outputs sit uncompressed as *.out outside any archive and are accepted by" % nplain,
          "   cs6_audit.read(). They are COUNTED AND NOT CLASSIFIED. Whether a file that is not a committed",
          "   archive member belongs in this corpus is a separate decision and is not taken here.",
          "",
          "(b) THE DISGUISED CHARGED RUNS - the reading that decides whether the audit was too permissive",
          "",
          "   A run whose charges are Ne> centres reports no point-charge file, so cs6_audit.py counted it as a",
          "   CHARGE-FREE CONTROL and let it set the thresholds every charged run was judged against."]
    dis = sorted([d for d in R if d["npc"] == 0 and d["nne"] > 0], key=lambda d: (d["arc"], d["run"]))
    if not dis:
        O += ["", "   There are none. E4.6's concern does not arise in this corpus, and the charge-detection fix",
              "   changes nothing. That is a result, not an absence of one."]
    else:
        O += ["",
              "   %d such runs, of which %d were inside archives and so WERE in D6's control pool."
              % (len(dis), sum(1 for d in dis if d["src"] == "archive")),
              "",
              "   %-34s %-26s %5s %4s %8s %8s %6s"
              % ("archive or file", "run", "charges", "nb", "sum|q|", "max|q|", "found")]
        for d in dis[:capped]:
            O.append("   %-34s %-26s %5d %4s %8s %8s %6s"
                     % (d["arc"][-34:], d["run"].split("/")[-1][:26], d["nne"],
                        d["nb"] if d["nb"] else "-",
                        "%.3f" % d["sumq"] if d["sumq"] is not None else "-",
                        "%.2f" % d["maxq"] if d["maxq"] is not None else "-", d["src"]))
        if len(dis) > capped:
            O.append("   ... and %d more" % (len(dis) - capped))
    O += ["",
          "(c) WHAT THE CONTROL MEDIANS MOVED TO",
          "",
          "   A control median that moves means every threshold anchored to it moved with it.",
          "",
          "   %5s %-12s %6s %6s %9s %9s %10s %9s"
          % ("nb", "basis", "stage", "n", "med cyc", "med gap", "med maxLw", "med dip")]
    allnb = sorted({nb for S in (S1, S2, S3) for nb in S["CTL"]}, key=lambda x: (x is None, x))
    for nb in allnb:
        bs = [d["bas"] for d in S3["cf"] if d["nb"] == nb]
        bs = max(set(bs), key=bs.count) if bs else "?"
        for S, lab in ((S1, "1 (D6)"), (S2, "2"), (S3, "3")):
            c = S["CTL"].get(nb)
            if not c:
                continue
            O.append("   %5s %-12s %6s %6d %9.1f %9.3f %10.3f %9.2f"
                     % (nb if nb else "-", bs, lab, c["n"], c.get("cyc", 0), c.get("gap", 0),
                        c.get("maxlw", 0), c.get("dip", 0)))
    O += ["",
          "   SPECIFICITY CHECK at each stage - charge-free controls the classifier itself flags. Zero is the",
          "   only value that licenses the rest of the table, and it is reported whatever it is.",
          ""]
    for S, lab in ((S1, "1 (D6)"), (S2, "2"), (S3, "3")):
        nbad = sum(1 for d in S["cf"] if d["v"] != "CLEAN")
        O.append("   stage %-7s %4d of %d" % (lab, nbad, len(S["cf"])))
    O += ["",
          "(d) VERDICTS ON POINT-CHARGE SINGLE POINTS, BY BASIS, AT EACH STAGE",
          "",
          "   %5s %-12s %6s %7s %9s %8s %9s %11s"
          % ("nb", "basis", "stage", "runs", "clean", "polarised", "strained", "collapsed")]
    for nb in sorted({d["nb"] for d in S3["sp"]}, key=lambda x: (x is None, x)):
        for S, lab in ((S1, "1"), (S2, "2"), (S3, "3")):
            g = [d for d in S["sp"] if d["nb"] == nb]
            if not g:
                continue
            bs = [d["bas"] for d in g]
            O.append("   %5s %-12s %6s %7d %9d %8d %9d %11d"
                     % (nb if nb else "-", max(set(bs), key=bs.count), lab, len(g),
                        vcount(g, "CLEAN"), vcount(g, "POLARISED"),
                        vcount(g, "STRAINED"), vcount(g, "COLLAPSED")))
    O += ["",
          "(e) EVERY VERDICT THAT CHANGED, ON THE RUNS D6 ALREADY HAD",
          "",
          "   Restricted to runs present in stage 1, so that a changed verdict is a consequence of a corrected",
          "   threshold and not of a newly opened file."]
    m1 = {key(d): d for d in S1["good"]}
    ch = []
    for d in S3["good"]:
        a = m1.get(key(d))
        if a and a["v"] != d["v"]:
            ch.append((a, d))
    sp3 = {key(d) for d in S3["sp"]}
    if not ch:
        O += ["", "   None. The corrected controls and the corrected charge detection change no verdict that D6",
              "   printed. D6's classification of the runs it saw stands."]
    else:
        des = [t for t in ch if key(t[1]) in sp3]
        oth = [t for t in ch if key(t[1]) not in sp3]
        O += ["",
              "   %d verdicts change: %d on point-charge single points, which is the table every Stage 6 number"
              % (len(ch), len(des)),
              "   is read off, and %d on runs that are not. The second group is expected and is not a new" % len(oth),
              "   failure: a charge-free optimisation that stage 1 let into the control pool was judged against a",
              "   median it had itself set, and once it is taken out of that pool it is judged against the single",
              "   points instead. It is listed so that the first group cannot be inflated by it.", ""]
        for tag, grp in (("POINT-CHARGE SINGLE POINTS", des), ("EVERYTHING ELSE", oth)):
            if not grp:
                O += ["   %s: none." % tag, ""]
                continue
            O += ["   %s" % tag, "",
                  "   %-34s %-26s %5s %4s   %-10s -> %-10s  %s"
                  % ("archive or file", "run", "chg", "nb", "stage 1", "stage 3", "why now")]
            for a, d in sorted(grp, key=lambda t: (t[0]["v"], t[1]["v"], t[0]["arc"]))[:capped]:
                O.append("   %-34s %-26s %5d %4s   %-10s -> %-10s  %s"
                         % (d["arc"][-34:], d["run"].split("/")[-1][:26], d["pc"],
                            d["nb"] if d["nb"] else "-", a["v"], d["v"],
                            "; ".join(d["w"]) or "no indicator departs"))
            if len(grp) > capped:
                O.append("   ... and %d more" % (len(grp) - capped))
            O.append("")
    O += ["",
          "(f) AT HOW MANY CHARGES DOES def2-SVP FIRST GO WRONG",
          "",
          "   def2-SVP (nb 264) is the basis every committed Stage 6 design was evaluated at. D6 read the onset",
          "   off an incomplete corpus; this is the same reading off the corrected one.",
          "",
          "   %6s %8s %7s %10s %10s %11s %10s %10s"
          % ("stage", "charges", "runs", "med cyc", "med gap", "worst gap", "strained", "collapsed")]
    for S, lab in ((S1, "1"), (S2, "2"), (S3, "3")):
        for lo, hi, tag in ((1, 1, "1"), (2, 2, "2"), (3, 3, "3"), (4, 4, "4"), (5, 9, "5-9"),
                            (10, 13, "10-13"), (14, 29, "14-29"), (30, 10 ** 6, ">=30")):
            g = [d for d in S["sp"] if d["nb"] == 264 and lo <= d["pc"] <= hi]
            if not g:
                continue
            md = lambda k: (lambda v: "%.3f" % st.median(v) if v else "-")([x[k] for x in g if x[k] is not None])
            O.append("   %6s %8s %7d %10s %10s %11s %10d %10d"
                     % (lab, tag, len(g), md("cyc"), md("gap"),
                        "%.3f" % min([x["gap"] for x in g if x["gap"] is not None], default=0),
                        vcount(g, "STRAINED"), vcount(g, "COLLAPSED")))
    low = [d["pc"] for d in S3["sp"] if d["nb"] == 264 and d["v"] in ("STRAINED", "COLLAPSED")]
    O += ["",
          "   LOWEST CHARGE COUNT AT WHICH A def2-SVP SINGLE POINT IS STRAINED OR COLLAPSED, stage 3: %s"
          % (min(low) if low else "none in the corpus")]
    O += ["",
          "(g) THE CHARGE VECTOR OF EVERY RUN WHOSE CHARGES ARE RECOVERABLE",
          "",
          "   Two sources, both committed, neither inferred: Ne> centres printed in the run's own coordinate",
          "   block, and a point-charge .pc file committed beside the output. ORCA does not echo a .pc file into",
          "   its output, so a run whose .pc file was not committed has no recoverable magnitudes and none is",
          "   guessed for it. This is the table section 6 of STAGE6_ASSESSMENT2_20261026.md needs, and it is the",
          "   first time it has had more than eight rows.",
          ""]
    rec = [d for d in S3["sp"] if d["chg"]]
    norec = [d for d in S3["sp"] if not d["chg"]]
    O += ["   %d of %d point-charge single points have recoverable charges (%d from Ne> centres, %d from a"
          % (len(rec), len(S3["sp"]), sum(1 for d in rec if d["qsrc"] == "Ne>"),
             sum(1 for d in rec if d["qsrc"] == ".pc")),
          "   committed .pc file); %d do not, and are excluded from reading (h) for that reason alone." % len(norec),
          ""]
    # one row per distinct charge vector AND geometry, because that is what (h) actually tests
    sig = {}
    for d in rec:
        sig.setdefault((d["qid"], d["geo"], d["nb"]), []).append(d)
    mixed = [(k, g) for k, g in sig.items()
             if {("sick" if x["v"] in ("STRAINED", "COLLAPSED") else "healthy") for x in g} == {"sick", "healthy"}]
    if rec:
        O += ["   %d distinct (charge vector, geometry, basis) signatures over those %d runs. One row each,"
              % (len(sig), len(rec)),
              "   sorted by peak field. Nothing is truncated.",
              "",
              "   %4s %4s %8s %7s %8s %9s %8s %6s %-24s  %s"
              % ("n", "nb", "sum|q|", "max|q|", "sum q^2", "peak |E|", "r min", "runs", "verdicts", "example")]
        for k in sorted(sig, key=lambda k: -(sig[k][0]["Epk"] or 0)):
            g = sig[k]
            e = g[0]
            vs = {}
            for x in g:
                vs[x["v"]] = vs.get(x["v"], 0) + 1
            O.append("   %4d %4s %8.3f %7.2f %8.3f %9s %8s %6d %-24s  %s"
                     % (len(e["chg"]), e["nb"] if e["nb"] else "-", e["sumq"], e["maxq"], e["sq2"],
                        "%.3f" % e["Epk"] if e["Epk"] is not None else "-",
                        "%.3f" % e["rmin"] if e["rmin"] is not None else "-",
                        len(g), ", ".join("%s %d" % (a, b) for a, b in sorted(vs.items())),
                        e["run"].split("/")[-1][:26]))
        O += ["", "   peak |E|  the largest |E| at any real atom, in V/Angstrom, E = %.5f q / r^2 summed over" % KE,
              "             that run's charges at that run's own final coordinates. r min is the closest",
              "             charge-to-atom approach, in Angstrom."]
    O += ["",
          "(h) THE PRE-REGISTERED TEST OF THE PEAK-FIELD HYPOTHESIS",
          "",
          "   STAGE6_ERRATUM5_20261026.txt offers peak local field as the only scalar still standing after sum|q|,",
          "   max|q|, n and sum q^2 all failed to separate eight arrangements. D8.4(h) fixed the test and its",
          "   falsification criterion BEFORE this corpus was scanned:",
          "",
          "       healthy = CLEAN or POLARISED (D6.3 states in terms that POLARISED is not a fault)",
          "       sick    = STRAINED or COLLAPSED",
          "       the hypothesis SURVIVES only if max(peak |E| over healthy) < min(peak |E| over sick),",
          "       that is, only if a separating threshold exists. ANY overlap falsifies it.",
          "",
          "   The test is run per basis set, because a field threshold that needed a different value at each",
          "   basis would not be a property of the design.",
          ""]
    ran = False
    for nb in sorted({d["nb"] for d in rec}, key=lambda x: (x is None, x)):
        g = [d for d in rec if d["nb"] == nb]
        hh = [d["Epk"] for d in g if d["v"] in ("CLEAN", "POLARISED") and d["Epk"] is not None]
        ss = [d["Epk"] for d in g if d["v"] in ("STRAINED", "COLLAPSED") and d["Epk"] is not None]
        if not hh or not ss:
            O += ["   nb %-6s %d healthy and %d sick runs carry a field. One side is empty, so the test cannot"
                  % (nb if nb else "-", len(hh), len(ss)),
                  "            run at this basis and the hypothesis is neither supported nor falsified here.", ""]
            continue
        ran = True
        ov = sum(1 for x in hh if x >= min(ss))
        O += ["   nb %-6s healthy %4d runs, peak |E| %.3f to %.3f" % (nb if nb else "-", len(hh), min(hh), max(hh)),
              "            sick    %4d runs, peak |E| %.3f to %.3f" % (len(ss), min(ss), max(ss)),
              "            %s %d healthy runs sit at or above the lowest sick run."
              % ("SURVIVES AT THIS BASIS:" if max(hh) < min(ss) else "FALSIFIED AT THIS BASIS:", ov)]
        cc = [d["Epk"] for d in g if d["v"] == "CLEAN" and d["Epk"] is not None]
        if cc:
            O.append("            secondary, CLEAN only against sick: %.3f to %.3f against %.3f up - %s"
                     % (min(cc), max(cc), min(ss), "separates" if max(cc) < min(ss) else "DOES NOT SEPARATE"))
        O.append("")
    if not ran:
        O += ["   THE TEST DID NOT RUN AT ANY BASIS. That is not a result about the hypothesis; it is a result",
              "   about the record. It means no basis set in this corpus has both a healthy and a sick run whose",
              "   charges are recoverable, so the committed outputs cannot decide the question either way.", ""]
    O += ["   A SIGNATURE CARRYING BOTH VERDICTS WOULD SETTLE MORE THAN THE PEAK FIELD. A signature here is the",
          "   exact charge vector and the exact final geometry to four decimals, matched by sha256 at one basis",
          "   - not a scalar summary of them. If two runs agree on all three and are classified differently then",
          "   NO function of the charges and the geometry can separate the verdicts, because any such function",
          "   returns one value for both runs. The difference would have to come from something that is not the",
          "   design: the SCF path, the initial guess, the ORCA build. That would mean the classifier is partly",
          "   reading run-to-run variation, and it would have to be settled before any threshold is fitted to it.",
          ""]
    if mixed:
        O += ["   %d signature%s carr%s both verdicts."
              % (len(mixed), "" if len(mixed) == 1 else "s", "ies" if len(mixed) == 1 else "y"),
              "",
              "   %4s %8s %9s %6s  %s" % ("n", "sum|q|", "peak |E|", "runs", "verdicts and an example of each")]
        for k, g in sorted(mixed, key=lambda t: -(t[1][0]["Epk"] or 0)):
            e = g[0]
            ex = {}
            for x in g:
                ex.setdefault(x["v"], x["run"].split("/")[-1][:24])
            O.append("   %4d %8.3f %9s %6d  %s"
                     % (len(e["chg"]), e["sumq"], "%.3f" % e["Epk"] if e["Epk"] is not None else "-",
                        len(g), "; ".join("%s %s" % (a, b) for a, b in sorted(ex.items()))))
        O.append("")
    else:
        O += ["   There are none: no signature in this corpus carries both verdicts. The question is left to the",
              "   per-basis test above.", ""]
    O += ["",
          "(i) THE GEOMETRY OPTIMISATIONS AND THE QM/MM RUNS, REPORTED APART", ""]
    for S, lab in ((S2, "2"), (S3, "3")):
        for nm, g in (("charge-free optimisations", S["cfo"]), ("point-charge optimisations", S["go"]),
                      ("QM/MM", S["mm"])):
            if g:
                O.append("   stage %s  %-28s %4d   clean %d, polarised %d, strained %d, collapsed %d, no control %d"
                         % (lab, nm, len(g), vcount(g, "CLEAN"), vcount(g, "POLARISED"),
                            vcount(g, "STRAINED"), vcount(g, "COLLAPSED"), vcount(g, "NO CONTROL")))
    O += ["",
          "WHAT THIS DOES NOT TEST",
          "",
          "   Nothing was computed. This cannot say what the correct barrier is for any run, only which numbers",
          "   should not be quoted. It inherits D6.0's disclosure in full: the four thresholds were chosen AFTER",
          "   the Stage 6 outputs had been inspected, so this is not a blind test either, and D8 deliberately",
          "   does not re-derive them - changing the corpus and the rule in one step would make the comparison",
          "   uninterpretable. Whether those two unanchored criteria should be replaced by SCF stability analysis",
          "   and FOD is the next question and is not this one.",
          "",
          "   A CLEAN verdict is necessary and not sufficient. Erratum 4 E4.6 records that B1 added a Pauli wall",
          "   and the gap recovered from 0.550 to 3.637 eV while b stayed 32.6 per cent from its target against a",
          "   10 per cent criterion, so a healthy wavefunction does not make a design right.",
          "",
          "   Reading (h) tests one scalar on one corpus at the bases that corpus happens to use. A scalar that",
          "   survives here is a candidate for a constraint, not a constraint: it would still have to be tested",
          "   on arrangements chosen to break it, which do not exist yet. Two electrostatic predictors have",
          "   already been proposed and failed in this project.",
          "",
          "   The charge magnitudes in (g) are recoverable only for runs that carry their charges as Ne> centres.",
          "   For point-charge-file runs the output does not echo the magnitudes and none is inferred.", ""]
    return "\n".join(O) + "\n"


def main(res, code, out):
    res, code, out = Path(res), Path(code), Path(out)
    try:
        if out.resolve().is_relative_to(res.resolve()):
            stop("refusing to write inside the results repo: %s" % out)
    except AttributeError:                                       # python < 3.9
        if str(out.resolve()).startswith(str(res.resolve())):
            stop("refusing to write inside the results repo: %s" % out)
    cs6 = code / "phase2.2" / "cs6_audit.py"
    if not cs6.is_file():
        stop("cs6_audit.py not found at %s" % cs6)
    hits = [p for p in res.rglob("DIAGNOSTIC6_RESULT.txt") if ".git" not in p.parts]
    if len(hits) != 1:
        stop("expected exactly one committed DIAGNOSTIC6_RESULT.txt in %s, found %d" % (res, len(hits)))
    d6 = hits[0]

    print("   reading the corpus ...", flush=True)
    R, nloose, nplain = collect(res)
    if not R:
        stop("found no completed ORCA single points on this dianion anywhere in %s" % res)
    S1 = stage(R, 1)

    print("   gate: rebuilding DIAGNOSTIC6_RESULT.txt from stage 1 ...", flush=True)
    got, want = render_d6(S1), d6.read_text()
    if got != want:
        g, w = got.splitlines(), want.splitlines()
        for i in range(max(len(g), len(w))):
            a = g[i] if i < len(g) else "<end of rebuild>"
            b = w[i] if i < len(w) else "<end of committed file>"
            if a != b:
                print("STOP: the rebuild of DIAGNOSTIC6_RESULT.txt does not match the committed file.",
                      file=sys.stderr)
                print("      first difference at line %d of %s" % (i + 1, d6), file=sys.stderr)
                print("      committed: %s" % b, file=sys.stderr)
                print("      rebuilt  : %s" % a, file=sys.stderr)
                print("", file=sys.stderr)
                print("      THAT IS THE FINDING. Either the committed corpus is no longer the corpus D6 was",
                      file=sys.stderr)
                print("      run against, or this script's reproduction of D6.3 is not faithful. No later",
                      file=sys.stderr)
                print("      reading would be interpretable, so none is written.", file=sys.stderr)
                raise SystemExit(1)
        stop("the rebuild differs from the committed file in trailing whitespace only")
    print("   gate passed: the rebuild is byte-identical to the committed file.", flush=True)

    S2, S3 = stage(R, 2), stage(R, 3)
    txt = render_d8(S1, S2, S3, R, nloose, nplain, str(d6.relative_to(res)), sha(d6), sha(cs6))
    out.mkdir(parents=True, exist_ok=True)
    (out / "DIAGNOSTIC8_RESULT.txt").write_text(txt)
    print(txt)
    print("   wrote %s" % (out / "DIAGNOSTIC8_RESULT.txt"))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        stop("usage: cs8_audit.py <results-repo> <code-repo> <out-dir>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
