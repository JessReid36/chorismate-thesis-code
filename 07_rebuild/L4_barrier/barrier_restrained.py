#!/usr/bin/env python3
"""Layer 4 -- BARRIER under a fixed field WITH one-sided flat-bottom geometric restraints (NEB-TS).

Route B, Design 1. Identical to the validated barrier.py EXCEPT it adds two one-sided flat-bottom
harmonic walls to the OpenMM system, alongside the existing LJ wall, in the SAME force-injection block.
  - C1-C6 (0,12) forming: field pries OPEN in reactant (3.12->4.38). Wall inert r<=R0_C1C6, harmonic above.
  - O3-C4 (7,8) breaking: field over-stretches in product (2.90->4.83). Wall inert r<=R0_O3C4, harmonic above.
Walls at BARE endpoint lengths (block field damage only; no TS compression). FIELD_OFF=1 -> no-field control.
"""
import os, sys
import numpy as np

ORCADIR = "/home/apps2/ORCA/6.0.1"
SIMPLE  = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF"
BLOCKS  = "%cpcm epsilon 4.0 end\n%scf MaxIter 300 Guess PModel end"
NCORES  = 8
LJ = {"C": (3.40, 0.086), "N": (3.25, 0.170), "O": (2.96, 0.210), "H": (1.06, 0.016)}
HARTREE_KCAL = 627.5094740631
REACT_IDX = list(range(24))
HESS_ATOMS = [0, 7, 8, 12, 3, 17]

R0_C1C6 = 3.15
R0_O3C4 = 2.95
KWALL   = 50.0
FIELD_OFF = os.environ.get("FIELD_OFF", "0") == "1"


def build_qmmm(coords, mmcharges, els):
    import openmm
    from ash import Fragment, ORCATheory, OpenMMTheory, QMMMTheory
    frag = Fragment(coordsstring=coords, charge=-2, mult=1, conncalc=False)
    orca = ORCATheory(orcadir=ORCADIR, orcasimpleinput=SIMPLE, orcablocks=BLOCKS,
                      numcores=NCORES, autostart=False)
    mm = OpenMMTheory(fragment=frag, dummysystem=True, platform="CPU", numcores=1,
                      autoconstraints=None, rigidwater=False)
    ang = openmm.unit.angstrom; kcal = openmm.unit.kilocalorie_per_mole
    nm = openmm.unit.nanometer; kj = openmm.unit.kilojoule_per_mole

    lj = openmm.CustomBondForce("4*eps*((sig/r)^12 - (sig/r)^6)")
    lj.addPerBondParameter("sig"); lj.addPerBondParameter("eps"); lj.setUsesPeriodicBoundaryConditions(False)
    for i in range(24):
        si, ei = LJ[els[i]]
        for j in range(24, frag.numatoms):
            sj, ej = LJ[els[j]]
            sig = ((si + sj) / 2 * ang).value_in_unit(nm)
            eps = (((ei * ej) ** 0.5) * kcal).value_in_unit(kj)
            lj.addBond(i, j, [sig, eps])
    lj.setForceGroup(11); mm.system.addForce(lj)

    wall = openmm.CustomBondForce("0.5*kw*(max(0, r - r0))^2")
    wall.addPerBondParameter("kw"); wall.addPerBondParameter("r0")
    wall.setUsesPeriodicBoundaryConditions(False)
    kw_kj = ((KWALL * kcal).value_in_unit(kj)) / ((1.0 * ang).value_in_unit(nm))**2
    wall.addBond(0, 12, [kw_kj, (R0_C1C6 * ang).value_in_unit(nm)])
    wall.addBond(7,  8, [kw_kj, (R0_O3C4 * ang).value_in_unit(nm)])
    wall.setForceGroup(12); mm.system.addForce(wall)

    try:
        mm.simulation.context.reinitialize(preserveState=True)
    except Exception:
        pass

    charges = [0.0] * len(mmcharges) if FIELD_OFF else mmcharges
    qmmm = QMMMTheory(qm_theory=orca, mm_theory=mm, fragment=frag, qmatoms=REACT_IDX,
                      charges=charges, embedding="elstat", qm_charge=-2, qm_mult=1, numcores=NCORES)
    return frag, qmmm, mm


def coords_from(coordsfile):
    L = open(coordsfile).read().splitlines()
    n = int(L[0].split()[0]); body = L[2:2 + n]
    els = [ln.split()[0] for ln in body]
    return "\n".join(ln.rstrip() for ln in body), els, n


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: python3 barrier_restrained.py <design_name>   (env FIELD_OFF=1 for control)")
    name = sys.argv[1]
    from ash import Fragment, NEBTS
    coords_r, els, n = coords_from("%s_coords.xyz" % name)
    mmcharges = [float(x) for x in open("%s_charges.txt" % name).read().split(",")]
    frag_r, qmmm_r, mm = build_qmmm(coords_r, mmcharges, els)

    e_mm = mm.run(current_coords=frag_r.coords, elems=frag_r.elems, Grad=False)
    print(">> MM energy (LJ + walls) live-check: %s   [FIELD_OFF=%s]" % (str(e_mm), FIELD_OFF))

    react = Fragment(xyzfile="%s_reactant_relaxed.xyz" % name)
    prod  = Fragment(xyzfile="%s_product_relaxed.xyz"  % name)
    res = NEBTS(reactant=react, product=prod, theory=qmmm_r,
                images=10, CI=True,
                ActiveRegion=True, actatoms=REACT_IDX,
                interpolation="IDPP",
                maxiter=300, OptTS_maxiter=300,
                hessian_for_TS="partial", partial_hessian_atoms=HESS_ATOMS,
                runmode='serial', printlevel=1)
    saddle = (getattr(res, "saddlepoint_fragment", None)
              or getattr(res, "Saddlepoint_fragment", None)
              or getattr(res, "saddlepoint", None))
    if saddle is None:
        raise RuntimeError("NEBTS result has no saddle attribute; dir=%s" % [a for a in dir(res) if not a.startswith("_")])
    E_ts = float(saddle.energy)
    E_r  = float(react.energy) if react.energy is not None else None
    if E_r is None:
        _, qmmm_re, _ = build_qmmm(coords_r, mmcharges, els)
        E_r = qmmm_re.run(current_coords=react.coords, elems=react.elems, Grad=False)
    barrier = (E_ts - E_r) * HARTREE_KCAL
    tag = "restrained_control" if FIELD_OFF else "restrained_field"
    out = ["design=%s" % name,
           "variant=%s" % tag,
           "walls=C1C6<=%.2f,O3C4<=%.2f,k=%.0f_kcal_per_A2" % (R0_C1C6, R0_O3C4, KWALL),
           "method=NEBTS_IDPP",
           "E_reactant_Eh=%.8f" % E_r,
           "E_saddle_Eh=%.8f" % E_ts,
           "barrier_kcal=%.2f" % barrier,
           "bare_reference_kcal=17.47",
           "delta_vs_bare_kcal=%.2f" % (barrier - 17.47)]
    open("%s_%s_barrier.txt" % (name, tag), "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
