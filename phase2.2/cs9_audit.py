#!/usr/bin/env python3
"""Stage 6 DIAGNOSTIC 9: the wavefunction-health audit with the charges actually detected.

NO NEW CALCULATION AND NO NETWORK. Diagnostic 6 is not edited, Diagnostic 8 is not edited, and
neither pre-registration is amended. This is a third audit with a third result file.

WHY IT EXISTS. Diagnostic 8 corrected D6's walker and categories and confirmed both defects, but its
own charge detector matched nothing in this corpus: it looked for Ne> centres on unprefixed lines
with four numbers, and ORCA prints the charge only in the input echo, prefixed '|  30> '. So 101
runs that carry charges entered D8's stage 3 as charge-free CONTROLS. D8's pre-registered
specificity check caught it - 49 of 469 controls flagged, against D6.4(a)'s rule that zero is the
only value that licenses the table - and STAGE6_DIAGNOSTIC8_ERRATUM1_20261026.txt withdraws D8's
stage 3 and its readings (b), (g) and (h). D8's stages 1 and 2 stand and are re-checked here.

WHAT THE CORPUS ACTUALLY CONTAINS, established by reading it rather than by assuming:

  - 1634 runs place BARE point charges through %pointcharges "<file>.pc". No Pauli wall.
  - 91 runs place charges as Ne> centres carrying a NewECP that replaces 10 core electrons. ORCA's
    own words for these, printed in the output, are "coreless ECP center with (optional) point
    charge". The ECP supplies the short-range repulsion a bare point charge has none of.
  - 10 runs place Ne> centres with NO ECP: bare charges written into the coordinate block.
  - 101 of those 101 Ne> runs report ORCA's point-charge count as zero, which is why D6 and D8 both
    filed them as controls.
  - 193 runs are QM/MM (! QMMM, %qmmm), whose charges are an MM environment. D6.1 excludes these by
    the proxy "more than 200 point charges"; they are detected here by the directive instead, which
    is what D6.1 says it means.

THE ECP/BARE DISTINCTION IS NOT BOOKKEEPING. Finding 4 identifies the absence of a Pauli wall as the
defect, and Erratum 4 E4.6 records B1 adding one. A charge with an ECP and a charge without are two
different physical models, and this diagnostic keeps them apart everywhere instead of pooling them.

cs6_audit.read() and all six of its thresholds are IMPORTED from cs6_audit.py and used unchanged.
D6.3's rule is applied unchanged. The variable is the corpus and how its charges are read.

usage: cs9_audit.py <results-repo> <code-repo> <out-dir>
"""
import sys, re, tarfile, lzma, math, hashlib, collections, statistics as st
from pathlib import Path

# importing cs6_audit would otherwise drop a __pycache__ directory into the code repo and leave it
# dirty, which is both a side effect on a repo this script only reads and a thing that would make
# the next run refuse. Set before the import, not after.
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cs6_audit as D6                      # read(), D_LW, R_DIP, GAP_ABS, R_GAP, R_CYC, MM_NPC

KE = 14.39964548          # e^2 / (4 pi eps0), in eV * Angstrom; E = KE*q/r^2 gives V/Angstrom
ANG = "CARTESIAN COORDINATES (ANGSTROEM)"

RE_ECHO = re.compile(r"^\|\s*\d+>\s?(.*)$")
RE_NEIN = re.compile(r"^\s*Ne>\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*(\S*)\s*$")
RE_PCD = re.compile(r'%\s*pointcharges\s+"([^"]+)"', re.I)
RE_QMMM = re.compile(r"(^|\s)!.*\bQMMM\b|^\s*%\s*qmmm\b", re.I | re.M)
RE_ATOM = re.compile(r"^\s*([A-Za-z][A-Za-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")
RE_NEXYZ = re.compile(r"^\s*Ne>\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")
# ORCA states the ECP assignment itself, per centre, in its own output. The input can declare the
# same ECP inline after the centre or in a %basis block; both produce these lines, and an audit that
# parsed the input syntax instead would call a %basis-declared run bare. It is what the program did
# that matters, not what the input looked like.
RE_ECPAT = re.compile(r"^Atom\s+\d+Ne>\s+ECP group\s*=>\s*\d+\s*$", re.M)
RE_LW = re.compile(r"^\s*\d+\s+([A-Za-z]{1,2})\s*:\s*(-?\d+\.\d+)")
# ORCA labels a point-charge centre "Q" in the Loewdin table and prints the centre's CHARGE there.
# cs6_audit.read() matches that label like any element, so for a run whose charges are coordinate-
# block centres its "largest Loewdin charge on a light atom" is the design's own charge read back,
# and D6.3's MOVED test then fires on any |q| >= 0.25 above the control with no reference to the
# electron density at all. D6.3 says "on a LIGHT ATOM"; a Q centre is not an atom. Excluding it
# implements that sentence rather than changing it - but both readings are reported, never one.
PC_LABEL = "Q"


def loewdin(s):
    """Returns (largest |Loewdin| over real atoms, number of point-charge centres in the table)."""
    k = s.rfind("LOEWDIN ATOMIC CHARGES")
    if k < 0:
        return None, 0
    at, nq = [], 0
    for l in s[k:].splitlines()[2:]:
        m = RE_LW.match(l)
        if not m:
            break
        if m.group(1) == PC_LABEL:
            nq += 1
        elif m.group(1) not in D6.METALS:
            at.append(abs(float(m.group(2))))
    return (max(at) if at else None), nq


def stop(m):
    print("STOP: %s" % m, file=sys.stderr)
    raise SystemExit(1)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def echo(s):
    """ORCA echoes the whole input file back, one line per '|  N> ' prefix. That echo is inside
    every output, so it is readable for an archived run exactly as for a loose one."""
    return "\n".join(m.group(1) for m in (RE_ECHO.match(l) for l in s.splitlines()) if m)


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


def geometry(s):
    """Real atoms, and Ne> centre positions, from the LAST Angstroem coordinate block. ORCA prints
    an Ne> centre there with its position only - the charge is in the echo, never here."""
    k = s.rfind(ANG)
    atoms, nepos = [], []
    if k < 0:
        return atoms, nepos
    for l in s[k:].splitlines()[2:]:
        if not l.strip():
            if atoms or nepos:
                break
            continue
        m = RE_NEXYZ.match(l)
        if m:
            nepos.append(tuple(float(g) for g in m.groups()))
            continue
        m = RE_ATOM.match(l)
        if m:
            atoms.append((float(m.group(2)), float(m.group(3)), float(m.group(4))))
            continue
        break
    return atoms, nepos


def field(atoms, chg):
    """Peak |E| over the real atoms, V/Angstrom, and the closest charge-to-atom approach."""
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
def collect(root):
    """One pass over everything, with each run's charges resolved to the best committed source."""
    root = Path(root)
    if not (root / ".git").is_dir():
        stop("%s is not a git repo" % root)

    disk_pc = collections.defaultdict(list)
    for p in root.rglob("*.pc"):
        if ".git" not in p.parts:
            disk_pc[p.name].append(p)

    R, nloose, nplain = [], 0, 0

    def take(s, where, run, src, arcpc=None):
        d = D6.read(s)
        if not d:
            return 0
        inp = echo(s)
        atoms, nepos = geometry(s)
        mlw_at, nq = loewdin(s)
        d.update(arc=where, run=run, src=src, nat=len(atoms), maxlw_at=mlw_at, nQ=nq,
                 qmmm=bool(RE_QMMM.search(inp)), note="")
        ne = []
        for l in inp.splitlines():
            m = RE_NEIN.match(l)
            if m:
                ne.append(tuple(float(m.group(i)) for i in (1, 2, 3, 4)))
        n_ecp_centres = len(RE_ECPAT.findall(s))
        if n_ecp_centres >= len(ne):
            ne_ecp, ne_bare = list(ne), []
        elif n_ecp_centres == 0:
            ne_ecp, ne_bare = [], list(ne)
        else:                      # some shielded and some not: say so rather than pick
            ne_ecp, ne_bare = list(ne), []
            d["note"] = "%d Ne> centres but ORCA assigns an ECP to only %d" % (len(ne), n_ecp_centres)
        chg, tier = [], "none"
        if ne:
            # the echo gives input coordinates; the coordinate block gives the final ones. For a
            # charge they must agree, and if they do not the field is not computed from a guess.
            ok = len(nepos) == len(ne)
            if ok:
                for (q, x, y, z) in ne:
                    if not any(abs(x - a) < 1e-3 and abs(y - b) < 1e-3 and abs(z - c) < 1e-3
                               for (a, b, c) in nepos):
                        ok = False
                        break
            if ok:
                chg, tier = ne, "echo"
            else:
                d["note"] = "Ne> echo and coordinate block disagree (%d vs %d centres)" % (len(ne), len(nepos))
        pcf = RE_PCD.search(inp)
        pcname = pcf.group(1).split("/")[-1] if pcf else None
        if pcname and not chg:
            text, how = None, None
            if src == "archive" and arcpc:
                dr = run.rsplit("/", 1)[0] if "/" in run else ""
                for k, v in arcpc.items():
                    if (k.rsplit("/", 1)[0] if "/" in k else "") == dr and k.rsplit("/", 1)[-1] == pcname:
                        text, how = v, "beside"
                        break
            elif src != "archive":
                cand = (root / where).parent / pcname
                if cand.is_file():
                    try:
                        text, how = cand.read_text(errors="replace"), "beside"
                    except Exception:
                        pass
            if text is None:
                hits = disk_pc.get(pcname, [])
                if len(hits) == 1:
                    try:
                        text, how = hits[0].read_text(errors="replace"), "named"
                    except Exception:
                        pass
            if text is not None:
                pc = read_pc(text)
                if pc and (d["npc"] == 0 or len(pc) == d["npc"]):
                    chg, tier = pc, how
                elif pc:
                    d["note"] = "%s holds %d charges, ORCA read %d" % (pcname, len(pc), d["npc"])
        d["ne_ecp"], d["ne_bare"], d["pcname"] = len(ne_ecp), len(ne_bare), pcname
        d["nchg"] = max(d["npc"], len(ne))
        d["nz"] = sum(1 for q, _, _, _ in ne if q != 0.0)
        if d["qmmm"]:
            d["cls"] = "qmmm"
        elif ne_ecp and d["nz"] == 0 and not pcname and not d["npc"]:
            # an ECP wall with every charge exactly zero. Not a control - something has been added
            # to the system - and not a charged design either. The wall on its own.
            d["cls"] = "wall"
        elif ne_ecp and (ne_bare or pcname or d["npc"]):
            d["cls"] = "mixed"
        elif ne_ecp:
            d["cls"] = "ecp"
        elif ne_bare or pcname or d["npc"]:
            d["cls"] = "bare"
        else:
            d["cls"] = "none"
        d["chg"], d["tier"] = chg, tier
        if chg:
            d["sumq"] = sum(abs(q) for q, _, _, _ in chg)
            d["maxq"] = max(abs(q) for q, _, _, _ in chg)
            d["sq2"] = sum(q * q for q, _, _, _ in chg)
            d["Epk"], d["rmin"] = field(atoms, chg)
            d["qid"] = hashlib.sha256(("|".join("%.4f,%.4f,%.4f,%.4f" % t
                                                for t in sorted(chg))).encode()).hexdigest()[:12]
        else:
            d["sumq"] = d["maxq"] = d["sq2"] = d["Epk"] = d["rmin"] = d["qid"] = None
        d["geo"] = hashlib.sha256(("|".join("%.4f,%.4f,%.4f" % t
                                            for t in atoms)).encode()).hexdigest()[:12] if atoms else "-"
        R.append(d)
        return 1

    for a in sorted(p for p in root.rglob("*.tar.xz") if ".git" not in p.parts):
        try:
            t = tarfile.open(str(a))
            members = [m for m in t.getmembers() if m.isfile()]
        except Exception:
            continue
        arcpc = {}
        for m in members:
            if m.name.endswith(".pc"):
                try:
                    arcpc[m.name] = t.extractfile(m).read().decode(errors="replace")
                except Exception:
                    pass
        for m in members:
            if not m.name.endswith(".out") or m.name.endswith(".pbs.out"):
                continue
            try:
                s = t.extractfile(m).read().decode(errors="replace")
            except Exception:
                continue
            take(s, str(a.relative_to(root)), m.name, "archive", arcpc)

    for a in sorted(p for p in root.rglob("*.xz")
                    if ".git" not in p.parts and not p.name.endswith(".tar.xz")):
        try:
            with lzma.open(str(a)) as fh:
                s = fh.read().decode(errors="replace")
        except Exception:
            continue
        if "FINAL SINGLE POINT ENERGY" not in s:
            continue
        nloose += take(s, str(a.relative_to(root)), a.name, "loose")

    for p in sorted(root.rglob("*.out")):
        if ".git" in p.parts or p.name.endswith(".pbs.out"):
            continue
        try:
            nplain += take(p.read_text(errors="replace"), str(p.relative_to(root)), p.name, "plain")
        except Exception:
            continue
    return R, nloose, nplain


# --------------------------------------------------------------------- D6.3, reproduced exactly
def verdict(d, CTL, atoms_only=False):
    """atoms_only applies D6.3 with the point-charge centres taken out of the Loewdin table, and
    with the dipole test suspended for runs that have them: ORCA's printed dipole for such a run
    includes those centres' own charges, so it is not the electron density's dipole and cannot be
    compared with a charge-free control. For every other run the two readings are identical."""
    c = CTL.get(d["nb"])
    if not c or "gap" not in c:
        return "NO CONTROL", []
    moved, sick = [], []
    lw = d["maxlw_at"] if atoms_only else d["maxlw"]
    if lw is not None and "maxlw" in c and lw - c["maxlw"] >= D6.D_LW:
        moved.append("Loewdin +%.3f e over control" % (lw - c["maxlw"]))
    skip_dip = atoms_only and d.get("nQ")
    if not skip_dip and d["dip"] is not None and c.get("dip") and d["dip"] / c["dip"] >= D6.R_DIP:
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
    """n = 1 D6 as published (gate 1) | 2 D8 stage 2 (gate 2) | 3 corrected."""
    pool = [dict(d) for d in R if d["src"] == "archive" or (n >= 3 and d["src"] == "loose")]
    for d in pool:
        d["pc"] = d["nchg"] if n >= 3 else d["npc"]
    good = [d for d in pool if d["norm"] and not d["notconv"]]
    CTL = {}
    for nb in {d["nb"] for d in good}:
        if n == 1:
            g = [d for d in good if d["nb"] == nb and d["pc"] == 0]
        elif n == 2:
            g = [d for d in good if d["nb"] == nb and d["pc"] == 0 and d["nopt"] == 0]
        else:
            g = [d for d in good if d["nb"] == nb and d["cls"] == "none" and d["nopt"] == 0]
        if g:
            CTL[nb] = {k: st.median([x[k] for x in g if x[k] is not None])
                       for k in ("cyc", "gap", "maxlw", "dip") if any(x[k] is not None for x in g)}
            CTL[nb]["n"] = len(g)
    for d in good:
        d["v"], d["w"] = verdict(d, CTL)
        d["v2"], d["w2"] = verdict(d, CTL, atoms_only=True) if n >= 3 else (d["v"], d["w"])
    if n == 1:
        cf = [d for d in good if d["pc"] == 0]
        sp = [d for d in good if 0 < d["pc"] <= D6.MM_NPC and d["nopt"] == 0]
        mm = [d for d in good if d["pc"] > D6.MM_NPC]
        go = [d for d in good if d["pc"] > 0 and d["nopt"] > 0]
        cfo, wall = [], []
    elif n == 2:
        cf = [d for d in good if d["pc"] == 0 and d["nopt"] == 0]
        cfo = [d for d in good if d["pc"] == 0 and d["nopt"] > 0]
        sp = [d for d in good if 0 < d["pc"] <= D6.MM_NPC and d["nopt"] == 0]
        go = [d for d in good if 0 < d["pc"] <= D6.MM_NPC and d["nopt"] > 0]
        mm = [d for d in good if d["pc"] > D6.MM_NPC]
        wall = []
    else:
        # QM/MM by its own directive, which is what D6.1 says the 200-charge proxy stands for
        mm = [d for d in good if d["cls"] == "qmmm"]
        rest = [d for d in good if d["cls"] != "qmmm"]
        cf = [d for d in rest if d["cls"] == "none" and d["nopt"] == 0]
        cfo = [d for d in rest if d["cls"] == "none" and d["nopt"] > 0]
        wall = [d for d in rest if d["cls"] == "wall"]
        sp = [d for d in rest if d["cls"] not in ("none", "wall") and d["nopt"] == 0]
        go = [d for d in rest if d["cls"] not in ("none", "wall") and d["nopt"] > 0]
    if n < 3:
        wall = []
    return dict(n=n, pool=pool, good=good, CTL=CTL, cf=cf, cfo=cfo, sp=sp, mm=mm, go=go, wall=wall)
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


# --------------------------------------------------------------------- gate 2: D8's licensed part
def d8_lines(S, lab_a, lab_d):
    """The exact lines D8 printed for this stage in its tables (a) and (d)."""
    out = ["   %-9s %7d %9d %8d %8d %8d %7d"
           % (lab_a, len(S["good"]), len(S["cf"]), len(S["cfo"]), len(S["sp"]),
              len(S["go"]), len(S["mm"]))]
    for nb in sorted({d["nb"] for d in S["sp"]}, key=lambda x: (x is None, x)):
        g = [d for d in S["sp"] if d["nb"] == nb]
        if not g:
            continue
        bs = [d["bas"] for d in g]
        out.append("   %5s %-12s %6s %7d %9d %8d %9d %11d"
                   % (nb if nb else "-", max(set(bs), key=bs.count), lab_d, len(g),
                      sum(1 for d in g if d["v"] == "CLEAN"), sum(1 for d in g if d["v"] == "POLARISED"),
                      sum(1 for d in g if d["v"] == "STRAINED"), sum(1 for d in g if d["v"] == "COLLAPSED")))
    return out


# --------------------------------------------------------------------- the readings
def vc(g, v):
    return sum(1 for d in g if d["v"] == v)


def render_d9(S1, S2, S3, R, nloose, nplain, d6p, d6s, d8p, d8s, cs6s, cap=120):
    O = ["STAGE 6 - DIAGNOSTIC 9 RESULT: the audit with the charges actually detected",
         "",
         "Measured by cs9_audit.py. NO NEW CALCULATION AND NO NETWORK. DIAGNOSTIC6_RESULT.txt and",
         "DIAGNOSTIC8_RESULT.txt are not edited and neither pre-registration is amended. D6.3's rule is",
         "applied unchanged, with every threshold imported from cs6_audit.py. What changes is the corpus and",
         "how its charges are read.",
         "",
         "   cs6_audit.py              sha256 %s" % cs6s,
         "   DIAGNOSTIC6_RESULT.txt    sha256 %s   (%s)" % (d6s, d6p),
         "   DIAGNOSTIC8_RESULT.txt    sha256 %s   (%s)" % (d8s, d8p),
         "",
         "THE TWO GATES, BOTH PASSED",
         "",
         "   GATE 1. The entire text of DIAGNOSTIC6_RESULT.txt was rebuilt from stage 1 and is byte-identical",
         "   to the committed file. Every number D6 printed comes back out of this script.",
         "   GATE 2. Stage 2 here reproduces, line for line, the stage 1 and stage 2 rows of DIAGNOSTIC8",
         "   tables (a) and (d) - the part of D8 that its own specificity check licenses and that",
         "   STAGE6_DIAGNOSTIC8_ERRATUM1_20261026.txt does NOT withdraw. D9 therefore agrees with everything",
         "   D8 is still entitled to say, and differs from it only where the erratum withdrew it.",
         "",
         "THE THREE STAGES",
         "",
         "   stage 1   D6 as published.                                                       gate 1",
         "   stage 2   D8 stage 2: controls restricted to charge-free single points, categories disjoint.",
         "                                                                                    gate 2",
         "   stage 3   corrected. The full corpus, charges read from the input echo and from committed .pc",
         "             files, QM/MM identified by its own directive, and ECP-shielded charges kept apart",
         "             from bare ones.",
         "",
         "(a) THE CORPUS, AND HOW ITS CHARGES ARE CARRIED",
         "",
         "   %-9s %7s %9s %8s %8s %8s %7s"
         % ("stage", "kept", "controls", "cf-opt", "chg-sp", "chg-opt", "QM/MM")]
    for S, lab in ((S1, "1 (D6)"), (S2, "2"), (S3, "3")):
        O.append("   %-9s %7d %9d %8d %8d %8d %7d"
                 % (lab, len(S["good"]), len(S["cf"]), len(S["cfo"]), len(S["sp"]),
                    len(S["go"]), len(S["mm"])))
    cls = collections.Counter(d["cls"] for d in R if d["norm"] and not d["notconv"])
    O += ["",
          "   %d completed runs on this dianion in all: %d inside archives, %d in loose .xz, %d as plain .out"
          % (sum(1 for d in R if d["norm"] and not d["notconv"]),
             sum(1 for d in R if d["norm"] and not d["notconv"] and d["src"] == "archive"),
             nloose, nplain),
          "   outside any archive. The plain .out runs are COUNTED AND NOT CLASSIFIED, as in D8: whether a",
          "   file that is not a committed archive member belongs in this corpus is a separate decision.",
          "",
          "   how the charges are carried, over every completed run:", ""]
    names = {"bare": "bare point charges, no Pauli wall (.pc file, or Ne> with no ECP)",
             "ecp": "ECP-shielded charges - a charge behind a Pauli wall",
             "wall": "the Pauli wall alone: ECP centres, every charge exactly zero",
             "mixed": "both kinds in one run",
             "qmmm": "QM/MM - an MM environment, not a design", "none": "nothing added at all"}
    for k, n in cls.most_common():
        O.append("   %-48s %6d" % (names.get(k, k), n))
    O += ["",
          "(b) THE RUNS D6 AND D8 FILED AS CONTROLS THAT ACTUALLY CARRY CHARGES",
          "",
          "   This is the reading D8 could not make. A charge placed as an Ne> centre is invisible to ORCA's",
          "   point-charge counter, so both earlier audits put these runs in the control pool, where they set",
          "   the thresholds every charged run was judged against."]
    bad = sorted([d for d in R if d["norm"] and not d["notconv"] and d["npc"] == 0
                  and (d["ne_ecp"] or d["ne_bare"])], key=lambda d: (d["arc"], d["run"]))
    if not bad:
        O += ["", "   There are none."]
    else:
        inD8 = sum(1 for d in bad if d["src"] in ("archive", "loose"))
        O += ["",
              "   %d such runs, %d of them inside the corpus D8's stage 3 classified." % (len(bad), inD8),
              "",
              "   %-40s %5s %5s %4s %10s %10s  %s"
              % ("file", "n", "kind", "nb", "sum|q|", "max|q|", "charges")]
        for d in bad[:cap]:
            O.append("   %-40s %5d %5s %4s %10s %10s  %s"
                     % (d["arc"][-40:], d["ne_ecp"] + d["ne_bare"],
                        "ECP" if d["ne_ecp"] else "bare", d["nb"] if d["nb"] else "-",
                        "%.5f" % d["sumq"] if d["sumq"] is not None else "-",
                        "%.5f" % d["maxq"] if d["maxq"] is not None else "-",
                        ", ".join("%+g" % q for q, _, _, _ in d["chg"][:6]) if d["chg"] else d["note"] or "-"))
        if len(bad) > cap:
            O.append("   ... and %d more" % (len(bad) - cap))
    O += ["",
          "(c) THE CONTROL MEDIANS, AND THE SPECIFICITY CHECK THAT LICENSES THEM",
          "",
          "   D6.4(a): zero is the only value that licenses the rest of the table. D8's stage 3 read 49 of 469",
          "   and is withdrawn for exactly that reason. This is the same check after the charges are detected.",
          "",
          "   %5s %-12s %6s %6s %9s %9s %10s %9s"
          % ("nb", "basis", "stage", "n", "med cyc", "med gap", "med maxLw", "med dip")]
    for nb in sorted({nb for S in (S1, S2, S3) for nb in S["CTL"]}, key=lambda x: (x is None, x)):
        bs = [d["bas"] for d in S3["cf"] if d["nb"] == nb]
        bs = max(set(bs), key=bs.count) if bs else "?"
        for S, lab in ((S1, "1 (D6)"), (S2, "2"), (S3, "3")):
            c = S["CTL"].get(nb)
            if c:
                O.append("   %5s %-12s %6s %6d %9.1f %9.3f %10.3f %9.2f"
                         % (nb if nb else "-", bs, lab, c["n"], c.get("cyc", 0), c.get("gap", 0),
                            c.get("maxlw", 0), c.get("dip", 0)))
    O += [""]
    lic = None
    for S, lab in ((S1, "1 (D6)"), (S2, "2"), (S3, "3")):
        nbad = sum(1 for d in S["cf"] if d["v"] != "CLEAN")
        O.append("   SPECIFICITY  stage %-7s %4d of %-5d  %s"
                 % (lab, nbad, len(S["cf"]),
                    "licensed" if nbad == 0 else "NOT LICENSED - the thresholds are too tight for this pool"))
        if S is S3:
            lic = (nbad == 0)
    if not lic:
        O += ["",
              "   STAGE 3 IS NOT LICENSED BY D6.4(a). Everything below that depends on stage 3's control",
              "   medians is reported for inspection and MUST NOT BE QUOTED as a classification. The readings",
              "   that do not depend on them - (a), (b), and the charge vectors in (g) - stand regardless."]
    O += ["",
          "(d) VERDICTS ON CHARGED SINGLE POINTS, BY BASIS AND BY MODEL CLASS, STAGE 3",
          "",
          "   Bare and ECP-shielded charges are two different physical models and are never pooled here.",
          "",
          "   TWO READINGS ARE GIVEN WHEREVER THEY DIFFER. ORCA labels a point-charge centre 'Q' in the",
          "   Loewdin table and prints the centre's own CHARGE on that line. cs6_audit.read() matches that",
          "   label like an element, so for a run whose charges sit in the coordinate block, D6.3's 'largest",
          "   Loewdin charge on a light atom' IS THE DESIGN'S OWN CHARGE READ BACK - and the MOVED test then",
          "   fires on any charge 0.25 e above the control with no reference to the electron density. The",
          "   second row takes those centres out, which is what D6.3's own words say, and suspends the dipole",
          "   test for those runs, because ORCA's printed dipole includes the centres' charges too. Runs whose",
          "   charges come from a .pc file are NOT affected: ORCA does not list them in the Loewdin table, so",
          "   the two readings coincide and every Stage 6 design verdict in D6 is untouched by this.",
          "",
          "   %5s %-12s %6s %7s %9s %8s %9s %11s  %s"
          % ("nb", "basis", "class", "runs", "clean", "polarised", "strained", "collapsed",
             "reading")]
    for nb in sorted({d["nb"] for d in S3["sp"]}, key=lambda x: (x is None, x)):
        for k in ("bare", "ecp", "mixed"):
            g = [d for d in S3["sp"] if d["nb"] == nb and d["cls"] == k]
            if not g:
                continue
            bs = [d["bas"] for d in g]
            nm = max(set(bs), key=bs.count)
            O.append("   %5s %-12s %6s %7d %9d %8d %9d %11d  %s"
                     % (nb if nb else "-", nm, k, len(g),
                        vc(g, "CLEAN"), vc(g, "POLARISED"), vc(g, "STRAINED"), vc(g, "COLLAPSED"),
                        "D6.3 as written"))
            if any(d["v"] != d["v2"] for d in g):
                O.append("   %5s %-12s %6s %7d %9d %8d %9d %11d  %s"
                         % ("", "", "", len(g),
                            sum(1 for d in g if d["v2"] == "CLEAN"),
                            sum(1 for d in g if d["v2"] == "POLARISED"),
                            sum(1 for d in g if d["v2"] == "STRAINED"),
                            sum(1 for d in g if d["v2"] == "COLLAPSED"),
                            "with the Q centres out of the Loewdin table (%d runs change)"
                            % sum(1 for d in g if d["v"] != d["v2"])))
    O += ["",
          "(e) EVERY VERDICT THAT CHANGED AGAINST D8 STAGE 2, ON THE RUNS BOTH HAD", ""]
    m2 = {(d["arc"], d["run"]): d for d in S2["good"]}
    ch = [(m2[(d["arc"], d["run"])], d) for d in S3["good"]
          if (d["arc"], d["run"]) in m2 and m2[(d["arc"], d["run"])]["v"] != d["v"]]
    if not ch:
        O += ["   None. Detecting the charges changes no verdict that D8 stage 2 printed."]
    else:
        O += ["   %d change." % len(ch), "",
              "   %-38s %-22s %5s %5s %4s   %-10s -> %-10s  %s"
              % ("file", "run", "n", "class", "nb", "D8 stage 2", "D9 stage 3", "why now")]
        for a, d in sorted(ch, key=lambda t: (t[0]["v"], t[1]["v"], t[0]["arc"]))[:cap]:
            O.append("   %-38s %-22s %5d %5s %4s   %-10s -> %-10s  %s"
                     % (d["arc"][-38:], d["run"].split("/")[-1][:22], d["pc"], d["cls"],
                        d["nb"] if d["nb"] else "-", a["v"], d["v"],
                        "; ".join(d["w"]) or "no indicator departs"))
        if len(ch) > cap:
            O.append("   ... and %d more" % (len(ch) - cap))
    O += ["",
          "(f) WHERE THE CHARGE VECTORS COME FROM",
          "",
          "   Three tiers, and the weakest is kept out of the test in (h).",
          "",
          "   echo     the charges are printed in the run's own input echo. Exact, no external file.",
          "   beside   a committed .pc file at the exact path the run's own directive names.",
          "   named    a committed .pc file matched only by a basename unique across the repository.",
          "            WEAKER. Reported, and excluded from (h).",
          "   none     the .pc file was never committed. No magnitude is guessed.",
          ""]
    for t in ("echo", "beside", "named", "none"):
        g = [d for d in S3["sp"] if d["tier"] == t]
        O.append("   %-8s %6d charged single points" % (t, len(g)))
    nt = [d for d in S3["sp"] if d["note"]]
    if nt:
        O += ["", "   %d runs were refused a charge vector by a consistency check rather than by absence:" % len(nt)]
        for d in nt[:12]:
            O.append("      %-44s %s" % (d["arc"][-44:], d["note"]))
    O += ["",
          "(g) THE CHARGE VECTORS, ONE ROW PER DISTINCT SIGNATURE",
          "",
          "   A signature is the exact charge vector and the exact final geometry to four decimals, matched",
          "   by sha256, at one basis and in one model class. Nothing is truncated.",
          ""]
    rec = [d for d in S3["sp"] if d["chg"] and d["tier"] in ("echo", "beside")]
    sig = {}
    for d in rec:
        sig.setdefault((d["qid"], d["geo"], d["nb"], d["cls"]), []).append(d)
    O += ["   %d charged single points with an exact charge vector, in %d distinct signatures."
          % (len(rec), len(sig)), ""]
    if rec:
        O.append("   %4s %5s %4s %8s %7s %8s %9s %8s %6s %-22s  %s"
                 % ("n", "class", "nb", "sum|q|", "max|q|", "sum q^2", "peak |E|", "r min",
                    "runs", "verdicts", "example"))
        for k in sorted(sig, key=lambda k: -(sig[k][0]["Epk"] or 0)):
            g = sig[k]; e = g[0]
            vs = collections.Counter(x["v"] for x in g)
            O.append("   %4d %5s %4s %8.3f %7.2f %8.3f %9s %8s %6d %-22s  %s"
                     % (len(e["chg"]), e["cls"], e["nb"] if e["nb"] else "-", e["sumq"], e["maxq"],
                        e["sq2"], "%.3f" % e["Epk"] if e["Epk"] is not None else "-",
                        "%.3f" % e["rmin"] if e["rmin"] is not None else "-", len(g),
                        ", ".join("%s %d" % (a, b) for a, b in sorted(vs.items())),
                        e["run"].split("/")[-1][:24]))
    O += ["",
          "(h) WHAT THE PAULI WALL ALONE DOES - the control series that was already run and never read",
          "",
          "   Finding 4 names the absence of a Pauli wall as the defect and Erratum 4 E4.6 records B1 adding",
          "   one. The corpus turns out to contain the matched series that question needs, in cs_stage4b and",
          "   cs_stage4c, and no diagnostic has read it:",
          "",
          "      nothing      a plain charge-free run. The control D6.3 anchors to.",
          "      wall only    Ne> ECP centres with every charge EXACTLY zero (the Z_* and *_q0 runs). The",
          "                   repulsive ECP is present and there is no charge. This isolates the wall.",
          "      wall+charge  the same ECP centres carrying charge.",
          "      bare charge  charge with no wall - what every Stage 6 design uses.",
          "",
          "   Each row is the median over runs of that class at that basis, against the same four indicators",
          "   D6.3 reads. The wall-only row is the one with no counterpart anywhere in the write-up.",
          "",
          "   %5s %-12s %-12s %6s %9s %9s %10s %9s %s"
          % ("nb", "basis", "class", "runs", "med cyc", "med gap", "med maxLw", "med dip", "verdicts")]
    pools = {"nothing": S3["cf"], "wall only": S3["wall"],
             "wall+charge": [d for d in S3["sp"] if d["cls"] == "ecp"],
             "bare charge": [d for d in S3["sp"] if d["cls"] == "bare"]}
    for nb in sorted({d["nb"] for g in pools.values() for d in g}, key=lambda x: (x is None, x)):
        shown = 0
        for lab in ("nothing", "wall only", "wall+charge", "bare charge"):
            g = [d for d in pools[lab] if d["nb"] == nb]
            if not g:
                continue
            shown += 1
            md = lambda k: (lambda v: "%.3f" % st.median(v) if v else "-")(
                [x[k] for x in g if x.get(k) is not None])
            vs = collections.Counter(d["v"] for d in g)
            bs = [d["bas"] for d in g]
            O.append("   %5s %-12s %-12s %6d %9s %9s %10s %9s %s"
                     % (nb if nb else "-", max(set(bs), key=bs.count), lab, len(g),
                        md("cyc"), md("gap"), md("maxlw"), md("dip"),
                        ", ".join("%s %d" % t for t in sorted(vs.items()))))
        if shown:
            O.append("")
    if not S3["wall"]:
        O += ["   There are no wall-only runs in this corpus, so the wall cannot be isolated from the charge",
              "   on this record.", ""]
    else:
        O += ["   READ THIS ROW WITH CARE. A wall-only run differs from a plain control by the ECP AND by the",
              "   three or five extra centres in the calculation. If its indicators sit on top of the control's,",
              "   the wall is inert on its own and any difference in the charged rows is the charge. If they do",
              "   not, the wall is not a free addition and B1's result cannot be read as the charge alone.",
              "   This is an observation on runs that already exist, not a designed experiment: the series were",
              "   run at different times for different purposes and are matched only in basis and geometry.", ""]

    O += ["",
          "(i) THE PEAK-FIELD HYPOTHESIS, TESTED WITHIN MODEL CLASS",
          "",
          "   D8.4(h)'s criterion, unchanged: healthy = CLEAN or POLARISED, sick = STRAINED or COLLAPSED, and",
          "   the hypothesis survives only if max(peak |E| over healthy) < min(peak |E| over sick). ANY overlap",
          "   falsifies it. D9 adds one restriction that D8 did not have and Erratum 5 did not observe: THE",
          "   TEST IS RUN WITHIN ONE MODEL CLASS. Comparing a bare charge with an ECP-shielded one varies the",
          "   Pauli wall and the field together, which is the confound",
          "   STAGE6_ERRATUM6_20261026.txt records in Erratum 5 E5.3.",
          ""]
    O += ["   The verdicts used here are the corrected ones from (d). The healthy/sick split is in fact",
          "   INVARIANT under that correction: taking the Q centres out of the Loewdin table can only move a",
          "   run between POLARISED and CLEAN, which are both healthy, or between COLLAPSED and STRAINED,",
          "   which are both sick. So this test gives the same answer on either reading.", ""]
    ran, verdicts = False, []
    for nb in sorted({d["nb"] for d in rec}, key=lambda x: (x is None, x)):
        for k in ("bare", "ecp", "mixed"):
            g = [d for d in rec if d["nb"] == nb and d["cls"] == k]
            hh = [d["Epk"] for d in g if d["v2"] in ("CLEAN", "POLARISED") and d["Epk"] is not None]
            ss = [d["Epk"] for d in g if d["v2"] in ("STRAINED", "COLLAPSED") and d["Epk"] is not None]
            if not g:
                continue
            if not hh or not ss:
                O.append("   nb %-5s %-6s %3d healthy, %3d sick - one side empty, the test does not run here."
                         % (nb if nb else "-", k, len(hh), len(ss)))
                continue
            ran = True
            ok = max(hh) < min(ss)
            verdicts.append((nb, k, ok))
            O += ["   nb %-5s %-6s healthy %4d runs, peak |E| %.3f to %.3f"
                  % (nb if nb else "-", k, len(hh), min(hh), max(hh)),
                  "                   sick    %4d runs, peak |E| %.3f to %.3f" % (len(ss), min(ss), max(ss)),
                  "                   %s  %d healthy runs sit at or above the lowest sick run."
                  % ("SURVIVES HERE:" if ok else "FALSIFIED HERE:",
                     sum(1 for x in hh if x >= min(ss)))]
    if verdicts:
        good_n = sum(1 for _, _, o in verdicts if o)
        O += ["",
              "   THE TEST RAN IN %d PLACES AND THE HYPOTHESIS SURVIVED IN %d OF THEM."
              % (len(verdicts), good_n)]
        if good_n == len(verdicts):
            O += ["   It is not falsified on this record. It is also not established: a scalar that separates",
                  "   where it has been looked at is a candidate for a constraint, not a constraint, and it has",
                  "   not been tested on arrangements chosen to break it."]
        else:
            O += ["   IT IS THEREFORE FALSIFIED. D8.4(h) fixed the criterion in advance and the test is run per",
                  "   basis precisely because a threshold needing a different value at each basis would not be a",
                  "   property of the design. Surviving at %d bases and failing at %d is that failure: there is no"
                  % (good_n, len(verdicts) - good_n),
                  "   single peak-field threshold that separates healthy from sick across this record.",
                  "   Peak local field joins sum|q|, max|q|, n and sum q^2. No scalar on the charge vector",
                  "   proposed so far survives contact with the corrected corpus."]
    if not ran:
        O += ["", "   THE TEST DID NOT RUN IN ANY CLASS AT ANY BASIS. That is a result about the record, not",
              "   about the hypothesis: no class and basis has both a healthy and a sick run with an exact",
              "   charge vector."]
    mixed = [(k, g) for k, g in sig.items()
             if len({"sick" if x["v"] in ("STRAINED", "COLLAPSED") else "healthy" for x in g}) == 2]
    O += ["",
          "   A signature carrying both verdicts would mean no function of the charges and the geometry can",
          "   separate them, and that the classifier is partly reading run-to-run variation.",
          "   %s" % ("%d signature(s) do." % len(mixed) if mixed else "There are none.")]
    for k, g in mixed[:20]:
        e = g[0]
        O.append("      n %d  %s  sum|q| %.3f  peak |E| %s  %s"
                 % (len(e["chg"]), e["cls"], e["sumq"],
                    "%.3f" % e["Epk"] if e["Epk"] is not None else "-",
                    ", ".join("%s %d" % t for t in sorted(collections.Counter(x["v"] for x in g).items()))))
    O += ["",
          "(j) THE OPTIMISATIONS AND THE QM/MM RUNS, REPORTED APART", ""]
    for nm, g in (("charge-free optimisations", S3["cfo"]), ("charged optimisations", S3["go"]),
                  ("QM/MM, by directive", S3["mm"])):
        if g:
            O.append("   %-28s %5d   clean %d, polarised %d, strained %d, collapsed %d, no control %d"
                     % (nm, len(g), vc(g, "CLEAN"), vc(g, "POLARISED"), vc(g, "STRAINED"),
                        vc(g, "COLLAPSED"), vc(g, "NO CONTROL")))
    O += ["",
          "WHAT THIS DOES NOT TEST",
          "",
          "   Nothing was computed. D6.0's disclosure is inherited unchanged: the four thresholds were chosen",
          "   after the Stage 6 outputs had been inspected, and D9 does not re-derive them, because changing",
          "   the corpus and the rule in one step would make the comparison uninterpretable. SCF stability",
          "   analysis and FOD remain the next question.",
          "",
          "   A CLEAN verdict is necessary and not sufficient. B1 recovered the gap from 0.550 to 3.637 eV with",
          "   a Pauli wall and b still stood 32.6 per cent from target against a 10 per cent criterion.",
          "",
          "   Reading (h) is observational. The bare and ECP-shielded families were run at different times for",
          "   different purposes and are matched only in basis and charge count. It can show that the wall",
          "   coincides with a different verdict distribution; it cannot show that the wall caused it.",
          "",
          "   Charge magnitudes exist only where the record carries them. For runs whose .pc file was never",
          "   committed there is no vector and none is inferred, and the size of that group is reported in (f)",
          "   because an untestable majority is itself a finding about the record.", ""]
    return "\n".join(O) + "\n"


def main(res, code, out):
    res, code, out = Path(res), Path(code), Path(out)
    if str(out.resolve()).startswith(str(res.resolve())):
        stop("refusing to write inside the results repo: %s" % out)
    cs6 = code / "phase2.2" / "cs6_audit.py"
    if not cs6.is_file():
        stop("cs6_audit.py not found at %s" % cs6)
    def one(name):
        h = [p for p in res.rglob(name) if ".git" not in p.parts]
        if len(h) != 1:
            stop("expected exactly one committed %s in %s, found %d" % (name, res, len(h)))
        return h[0]
    d6, d8 = one("DIAGNOSTIC6_RESULT.txt"), one("DIAGNOSTIC8_RESULT.txt")

    print("   reading the corpus ...", flush=True)
    R, nloose, nplain = collect(res)
    if not R:
        stop("found no completed ORCA single points on this dianion anywhere in %s" % res)
    S1 = stage(R, 1)

    print("   gate 1: rebuilding DIAGNOSTIC6_RESULT.txt from stage 1 ...", flush=True)
    got, want = render_d6(S1), d6.read_text()
    if got != want:
        g, w = got.splitlines(), want.splitlines()
        for i in range(max(len(g), len(w))):
            a = g[i] if i < len(g) else "<end of rebuild>"
            b = w[i] if i < len(w) else "<end of committed file>"
            if a != b:
                print("STOP: gate 1 failed at line %d of %s" % (i + 1, d6), file=sys.stderr)
                print("      committed: %s\n      rebuilt  : %s" % (b, a), file=sys.stderr)
                break
        stop("the rebuild of DIAGNOSTIC6_RESULT.txt does not match the committed file")
    print("   gate 1 passed.", flush=True)

    S2 = stage(R, 2)
    print("   gate 2: checking D8's licensed stage 1 and stage 2 rows ...", flush=True)
    d8txt = set(d8.read_text().splitlines())
    missing = [l for l in d8_lines(S1, "1 (D6)", "1") + d8_lines(S2, "2", "2") if l not in d8txt]
    if missing:
        print("STOP: gate 2 failed. These lines are not in the committed DIAGNOSTIC8_RESULT.txt:",
              file=sys.stderr)
        for l in missing[:10]:
            print("      %s" % l, file=sys.stderr)
        stop("D9 does not reproduce the part of D8 that its erratum leaves standing")
    print("   gate 2 passed.", flush=True)

    S3 = stage(R, 3)
    txt = render_d9(S1, S2, S3, R, nloose, nplain, str(d6.relative_to(res)), sha(d6),
                    str(d8.relative_to(res)), sha(d8), sha(cs6))
    out.mkdir(parents=True, exist_ok=True)
    (out / "DIAGNOSTIC9_RESULT.txt").write_text(txt)
    print(txt)
    print("   wrote %s" % (out / "DIAGNOSTIC9_RESULT.txt"))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        stop("usage: cs9_audit.py <results-repo> <code-repo> <out-dir>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
