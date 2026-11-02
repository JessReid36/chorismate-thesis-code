#!/usr/bin/env python3
"""cs6_exp0.py -- EXPERIMENT 0: the pseudopotential convention, and whether each
published parameter set reproduces its own fitting target.

Stage 6, 2 November 2026. No HPC time. Run as:

    python3 -IB cs6_exp0.py > EXPERIMENT0_RESULT.txt

Requires pyscf and basis_set_exchange (pip install pyscf basis_set_exchange).
Both ship their data locally; no network access is used.

WHAT THIS ANSWERS
-----------------
Marefat Khah et al. (acs.jctc.9b01162) print their pseudopotential parameters
"in TURBOMOLE input format" with 2 in the power column, and state that the
exponents n are set to zero. The project transcribes them into ORCA NewECP
blocks, whose 4th column the ORCA 6.0 manual calls "the radial power n_kl"
without giving an equation. If ORCA's column does not mean what TURBOMOLE's
means, every walled number in the Stage 6 record used a different potential
from the published one.

The SI's captions make the fit falsifiable: a centre carrying ZERO electrons
must reproduce named orbital energies of the neutral atom in a named basis. So
build T + V_ECP for such a centre, diagonalise it per angular momentum, and
compare with all-electron Hartree-Fock. Under the NWChem/PySCF convention the
power column p means r^(p-2), so p = 2 is r^0 and p = 4 is r^2. Run both.

READINGS
--------
  1  Does the Ne set reproduce its target under r^0 and not under r^2?
  2  Is the agreement specific to the doubly augmented basis the SI names?
  3  Do all three published sets reproduce their own targets?
  4  Is the ORCA power column the same column as NWChem's?
  5  Can Ar be validated on ANY available basis, under either reading of its
     d target, and does its coarser quoted precision explain the residual?
"""

import sys

import numpy as np

try:
    import basis_set_exchange as bse
    from pyscf import gto, scf
except ImportError as e:
    sys.stderr.write("STOP  %s\n      pip install pyscf basis_set_exchange\n" % e)
    sys.exit(2)

LMAP = "spdfghi"
H2KCAL = 627.5094740631


def say(s=""):
    sys.stdout.write(s + "\n")
    sys.stdout.flush()


def stop(code, msg):
    say("\nSTOP  %s\n      %s" % (code, msg))
    sys.exit(3)


# ---------------------------------------------------------------- SI Tables S1-S3
# Read from chorismate-refs/ct9b01162_si_001.pdf p.S2-S3. The SI prints
# coefficient / power / exponent; ORCA's order is exponent / coefficient / power.
# "target" is the shells the caption names; "digits" is the precision it claims.
SETS = {
    "He": dict(
        Z=2, ncore=2, lmax="d", basis="aug-cc-pV5Z", naug=1,
        target="2s2p3s3p", digits="within 6 digits",
        channels=[("S", 0.5098, 2.4200), ("P", 0.49165, -0.43590)],
        local=(1.0000, 0.0000),
        levels=[("s", 0), ("s", 1), ("p", 0), ("p", 1)],
    ),
    "Ne": dict(
        Z=10, ncore=10, lmax="f", basis="aug-cc-pwCV5Z", naug=1,
        target="3s3p3d4s4p4d", digits="within 6 digits",
        channels=[("S", 2.0475, 54.5100), ("P", 0.44815, 1.46500), ("D", 0.49205, -0.83800)],
        local=(1.0000, 0.0000),
        levels=[("s", 0), ("s", 1), ("p", 0), ("p", 1), ("d", 0), ("d", 1)],
    ),
    "Ar": dict(
        Z=18, ncore=18, lmax="f", basis="aug-cc-pV(5+d)Z", naug=1,
        target="4s4p4d5s5p5d", digits="NONE STATED",
        channels=[("S", 1.6410, 275.00), ("P", 0.2733, 1.9000), ("D", 0.44000, -3.4000)],
        local=(1.0000, 0.0000),
        levels=[("s", 0), ("s", 1), ("p", 0), ("p", 1), ("d", 0), ("d", 1)],
        note="the SI names aug-cc-pV(6+d)Z, which basis-set-exchange does not "
             "carry; aug-cc-pV(5+d)Z is the nearest available and is a substitute",
    ),
}

# ORCA 6.0 manual, section 7.5, its own worked example for V with SD(10,MDF).
# [quoted] -- read from the manual page, not computed here.
ORCA_V_EXAMPLE = [(1, "14.4900000000", "178.4479710000", 2),
                  (2, "6.5240000000", "19.8313750000", 2)]


def basis_for(el, name, naug):
    kw = {} if not naug else {"augment_diffuse": naug}
    return gto.basis.parse(
        bse.get_basis(name, elements=[el], fmt="nwchem", header=False, **kw))


def ecp_text(el, S, power):
    t = "%s nelec %d\n%s ul\n%d %.8f %.8f\n" % (el, S["ncore"], el, power,
                                                S["local"][0], S["local"][1])
    for ch, ex, co in S["channels"]:
        t += "%s %s\n%d %.8f %.8f\n" % (el, ch, power, ex, co)
    return t


def lblocks(mol):
    lb, off = {}, 0
    for b in range(mol.nbas):
        l = mol.bas_angular(b)
        n = mol.bas_nctr(b) * (2 * l + 1)
        lb.setdefault(LMAP[l], []).extend(range(off, off + n))
        off += n
    return {k: np.array(v) for k, v in lb.items()}


def assign_l(c, S, lb):
    """Dominant angular momentum of an eigenvector, by overlap-weighted weight."""
    best, bl = -1.0, None
    for l, idx in lb.items():
        w = float(c[idx] @ S[np.ix_(idx, idx)] @ c[idx])
        if w > best:
            best, bl = w, l
    return bl, best


def levels(vals, tol=1e-7):
    """Collapse a symmetry block into DISTINCT levels with their degeneracies.
    p has three degenerate components and d has five, so indexing the raw
    eigenvalue list returns the same level several times."""
    out = []
    for v in sorted(vals):
        if not out or abs(v - out[-1][0]) > tol:
            out.append([v, 1])
        else:
            out[-1][1] += 1
    return out


def reference(el, bas):
    """All-electron RHF; distinct unoccupied levels per angular momentum."""
    mol = gto.M(atom="%s 0 0 0" % el, basis={el: bas}, spin=0, verbose=0)
    mf = scf.RHF(mol).run()
    if not mf.converged:
        stop("EXP0.D1", "the all-electron RHF reference for %s did not converge" % el)
    S, lb = mol.intor("int1e_ovlp"), lblocks(mol)
    raw = {}
    for i in range(mol.nelectron // 2, mol.nao):
        l, w = assign_l(mf.mo_coeff[:, i], S, lb)
        if w >= 0.90:
            raw.setdefault(l, []).append(mf.mo_energy[i])
    return mol, mf, lb, {l: levels(v) for l, v in raw.items()}


def ecp_spectrum(el, bas, S, power, lb):
    """Distinct levels of T + V_ECP for a centre with no electrons and no charge."""
    em = gto.M(atom="%s 0 0 0" % el, basis={el: bas},
               ecp={el: gto.basis.parse_ecp(ecp_text(el, S, power))},
               spin=0, verbose=0)
    if em.nelectron != 0 or em.atom_charges()[0] != 0:
        stop("EXP0.G1", "the %s centre carries %d electron(s) and charge %s; the fit is "
                        "defined for a centre with neither."
             % (el, em.nelectron, em.atom_charges()))
    Vnuc = em.intor("int1e_nuc")
    if abs(Vnuc).max() > 1e-10:
        stop("EXP0.G1", "the %s centre's nuclear-attraction matrix is not zero "
                        "(max %.2e)" % (el, abs(Vnuc).max()))
    Sm = em.intor("int1e_ovlp")
    Hm = em.intor("int1e_kin") + em.intor("ECPscalar")
    Li = np.linalg.inv(np.linalg.cholesky(Sm))
    w, c = np.linalg.eigh(Li @ Hm @ Li.T)
    cao = Li.T @ c
    raw = {}
    for i in range(len(w)):
        l, wt = assign_l(cao[:, i], Sm, lb)
        if wt >= 0.90:
            raw.setdefault(l, []).append(w[i])
    return {l: levels(v) for l, v in raw.items()}


def deviations(ref, got, want):
    d = []
    for l, k in want:
        if l not in ref or len(ref[l]) <= k or l not in got or len(got[l]) <= k:
            stop("EXP0.D2", "level %s[%d] is not available on both sides" % (l, k))
        d.append((l, k, ref[l][k][0], got[l][k][0], got[l][k][0] - ref[l][k][0]))
    return d


def main():
    say("=" * 78)
    say("EXPERIMENT 0  --  the pseudopotential convention")
    say("Stage 6, 2 November 2026.  No HPC time.  pyscf %s" % __import__("pyscf").__version__)
    say("=" * 78)
    say("\nEvery number below is measured by this script. The ORCA manual values in")
    say("READING 4 are quoted from the manual and marked as such.")

    # ------------------------------------------------------ READING 1 and 3
    say("\n" + "=" * 78)
    say("READING 1 and 3: does each published set reproduce its own fitting target?")
    say("=" * 78)
    summary = {}
    for el in ("He", "Ne", "Ar"):
        S = SETS[el]
        say("\n--- %s  (SI target %s, basis %s, precision claimed: %s)"
            % (el, S["target"], S["basis"], S["digits"]))
        if "note" in S:
            say("    note: %s" % S["note"])
        bas = basis_for(el, S["basis"], S["naug"])
        mol, mf, lb, ref = reference(el, bas)
        say("    reference RHF  E = %.9f Eh,  nao %d" % (mf.e_tot, mol.nao))
        row = {}
        for power, label in ((2, "r^0"), (4, "r^2")):
            got = ecp_spectrum(el, bas, S, power, lb)
            dev = deviations(ref, got, S["levels"])
            row[label] = max(abs(x[4]) for x in dev)
            if power == 2:
                say("    %-4s level      HF reference        T + V_ECP        difference" % "")
                for l, k, r, g, d in dev:
                    say("         %s[%d]        %+.9f      %+.9f     %+.3e" % (l, k, r, g, d))
            say("    max |deviation| under %s : %.3e Eh" % (label, row[label]))
        summary[el] = row

    say("\n" + "-" * 78)
    say("   %-4s %-12s %-14s %-12s %-12s %s" %
        ("set", "target", "precision", "r^0", "r^2", "verdict"))
    for el in ("He", "Ne", "Ar"):
        r0, r2 = summary[el]["r^0"], summary[el]["r^2"]
        v = "VALIDATED" if r0 <= 1e-5 else "NOT VALIDATED"
        say("   %-4s %-12s %-14s %.3e    %.3e    %s"
            % (el, SETS[el]["target"], SETS[el]["digits"], r0, r2, v))

    # ------------------------------------------------------ READING 2
    say("\n" + "=" * 78)
    say("READING 2: is the Ne agreement specific to the basis the SI names?")
    say("=" * 78)
    S = SETS["Ne"]
    say("   %-8s %5s  %-12s %-12s" % ("aug level", "nao", "r^0", "r^2"))
    for naug, name in ((0, "aug"), (1, "d-aug"), (2, "t-aug")):
        bas = basis_for("Ne", S["basis"], naug)
        mol, mf, lb, ref = reference("Ne", bas)
        out = []
        for power in (2, 4):
            got = ecp_spectrum("Ne", bas, S, power, lb)
            out.append(max(abs(x[4]) for x in deviations(ref, got, S["levels"])))
        say("   %-8s %5d  %.3e    %.3e" % (name, mol.nao, out[0], out[1]))
    say("\n   The SI names d-aug-cc-pwCV5Z. If the agreement is sharply best there, the")
    say("   caption is meant as written and the even-tempered construction used here is")
    say("   close to the authors'.")

    # ------------------------------------------------------ READING 4
    say("\n" + "=" * 78)
    say("READING 4: is ORCA's radial-power column the same column as NWChem's?")
    say("=" * 78)
    say("   ORCA's own manual (6.0, section 7.5) illustrates NewECP with a published")
    say("   potential, V with SD(10,MDF). Its s channel reads [quoted]:")
    for ix, ex, co, p in ORCA_V_EXAMPLE:
        say("        %d  %s  %s  %d" % (ix, ex, co, p))
    # Collect EVERY listing that carries those primitives, not the first one found:
    # several independent sets agreeing is stronger than one, and the first hit
    # alphabetically may be a set that merely embeds the potential.
    found = []
    for name in bse.get_all_basis_names():
        try:
            txt = bse.get_basis(name, elements=["V"], fmt="nwchem", header=False)
        except Exception:
            continue
        flat = txt.replace(" ", "")
        if "14.49" in flat and "178.447971" in flat:
            found.append((name, txt))
    found.sort(key=lambda nt: (0 if "stuttgart" in nt[0].lower() else 1, nt[0]))
    hit = found[0] if found else None
    if hit is None:
        say("\n   basis-set-exchange carries no V ECP with those primitives, so this")
        say("   cross-check could not be made here. READING 4: NOT ESTABLISHED.")
    else:
        name, txt = hit
        say("\n   basis-set-exchange carries those primitives in %d set(s): %s"
            % (len(found), ", ".join(n for n, _ in found)))
        say("   Taking %s, which lists the potential in NWChem format:" % name)
        i = txt.upper().find("ECP")
        cols = set()
        for ln in txt[i:].splitlines():
            c = ln.split()
            if len(c) == 3 and (c[1].startswith("14.49") or c[1].startswith("6.524")):
                say("        %s  %s  %s" % (c[0], c[1], c[2]))
                cols.add(int(float(c[0])))
        orca_cols = {p for _, _, _, p in ORCA_V_EXAMPLE}
        if cols != orca_cols:
            stop("EXP0.G2", "the power column reads %s in the NWChem listing and %s in "
                            "ORCA's example; they are not the same column."
                 % (sorted(cols), sorted(orca_cols)))
        say("\n   NWChem documents its radial form as")
        say("        U_l(r) = r^(-2) * sum_k A_lk * r^(n_lk) * exp(-B_lk r^2)")
        say("   and states the shift is 'in accordance with most of the ECP literature")
        say("   and implementations', so n = 2 means r^0. [quoted, NWChem ECP manual]")
        say("\n   Same potential, same exponents, same coefficients, and the power column")
        say("   reads 2 in both listings. READING 4: ORCA's column is NWChem's column,")
        say("   so the project's blocks are r^0 blocks.")

    # ------------------------------------------------------ READING 5
    say("\n" + "=" * 78)
    say("READING 5: can Ar be validated on any available basis, and is its coarser")
    say("           quoted precision the explanation?")
    say("=" * 78)
    A = SETS["Ar"]
    SP = [("s", 0), ("s", 1), ("p", 0), ("p", 1)]
    say("   %-18s %-6s %5s  %-12s %-12s" % ("basis", "aug", "nao", "d = 3d,4d", "d = 4d,5d"))
    best = (1e9, None)
    for name in ("aug-cc-pV(5+d)Z", "aug-cc-pV6Z", "aug-cc-pV5Z", "aug-cc-pwCV5Z"):
        for naug, lab in ((0, "aug"), (1, "d-aug")):
            try:
                bas = basis_for("Ar", name, naug)
            except Exception:
                say("   %-18s %-6s  unavailable" % (name, lab))
                continue
            mol, mf, lb, ref = reference("Ar", bas)
            got = ecp_spectrum("Ar", bas, A, 2, lb)
            cells = []
            for dks in ((0, 1), (1, 2)):
                d = deviations(ref, got, SP + [("d", k) for k in dks])
                m = max(abs(x[4]) for x in d)
                cells.append(m)
                if m < best[0]:
                    best = (m, "%s %s, d=%s" % (name, lab, "3d,4d" if dks == (0, 1) else "4d,5d"))
            say("   %-18s %-6s %5d  %.3e    %.3e" % (name, lab, mol.nao, cells[0], cells[1]))
    say("\n   best cell: %.3e Eh at %s" % best)
    if best[0] <= 1e-5:
        say("   Ar CAN be validated on an available basis; the headline verdict is wrong.")
        stop("EXP0.G3", "READING 1 reported Ar as NOT VALIDATED but READING 5 found a basis "
                        "where it validates (%.3e at %s)." % best)
    say("   No available basis validates Ar. The SI names aug-cc-pV(6+d)Z, which")
    say("   basis-set-exchange does not carry.")

    say("\n   Is the coarser quoted precision the explanation? Ar's parameters carry 4-5")
    say("   significant figures where Ne's carry 6. Moving each coefficient and exponent")
    say("   within half a unit of its last quoted digit, 24 draws:")
    import random
    bas = basis_for("Ar", A["basis"], A["naug"])
    mol, mf, lb, ref = reference("Ar", bas)
    LSD = [("S", 0.00005, 0.005), ("P", 0.00005, 0.00005), ("D", 0.000005, 0.00005)]
    random.seed(7)
    lo, hi = 1e9, 0.0
    for _ in range(24):
        pert = dict(A)
        pert["channels"] = [(ch, ex + random.uniform(-dex, dex), co + random.uniform(-dco, dco))
                            for (ch, ex, co), (_, dex, dco) in zip(A["channels"], LSD)]
        got = ecp_spectrum("Ar", bas, pert, 2, lb)
        m = max(abs(x[4]) for x in deviations(ref, got, A["levels"]))
        lo, hi = min(lo, m), max(hi, m)
    say("      residual range %.3e .. %.3e Eh   (spread %.1e)" % (lo, hi, hi - lo))
    if hi - lo > 1e-4:
        say("      The spread is comparable to the residual, so truncation may explain it.")
    else:
        say("      The spread is %.0f times too small to account for a residual of %.3e,"
            % (lo / max(hi - lo, 1e-12), lo))
        say("      so TRUNCATED PRECISION IS EXCLUDED as the cause.")

    # ------------------------------------------------------ verdicts
    say("\n" + "=" * 78)
    say("VERDICTS")
    say("=" * 78)
    conv = summary["Ne"]["r^0"] <= 1e-5 and summary["Ne"]["r^2"] > 1e-2
    say("   CONVENTION %s" % ("CONFIRMED" if conv and hit else "NOT CONFIRMED"))
    for el in ("He", "Ne", "Ar"):
        say("   %-3s %s   (max |deviation| %.3e Eh under r^0)"
            % (el, "VALIDATED    " if summary[el]["r^0"] <= 1e-5 else "NOT VALIDATED",
               summary[el]["r^0"]))
    say("")
    if conv and hit:
        say("   The operator ORCA was given is the authors' operator. STAGE6_ERRATUM7")
        say("   E7.8's concern is about domain of applicability, not transcription.")
    if summary["Ar"]["r^0"] > 1e-5:
        say("   The Ar block cannot be validated against its own published target. Any")
        say("   experiment using it reports a sensitivity probe, not the published Ar")
        say("   operator. See ECP_CONVENTION_20261102.md section 3.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
