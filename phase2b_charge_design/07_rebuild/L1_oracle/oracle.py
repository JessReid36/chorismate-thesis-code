#!/usr/bin/env python3
"""Layer 1 -- the ORACLE (GOCAT-faithful evaluator). Clean-room rebuild."""
import os
import sys

ORCADIR = "/home/apps2/ORCA/6.0.1"
SIMPLE  = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF"
BLOCKS  = "%cpcm epsilon 4.0 end\n%scf MaxIter 300 Guess PModel end"
NCORES  = 8
LJ = {"C": (3.40, 0.086), "N": (3.25, 0.170), "O": (2.96, 0.210), "H": (1.06, 0.016)}
CONV = {"convergence_energy": 1e-5, "convergence_grms": 3.0e-3, "convergence_gmax": 6.0e-3,
        "convergence_drms": 1.0e-2, "convergence_dmax": 1.5e-2}
COV = {"C": 0.77, "O": 0.66, "N": 0.70, "H": 0.31}
BOND_TOL = 1.30
REACTING = {(7, 8), (0, 12)}
HARTREE_KCAL = 627.5094740631


def _read_body(path, a, b):
    return open(path).read().splitlines()[a:b]


def _bonds(els, R):
    s = set()
    for i in range(len(els)):
        for j in range(i + 1, len(els)):
            d = sum((R[i][k] - R[j][k]) ** 2 for k in range(3)) ** 0.5
            if d < BOND_TOL * (COV.get(els[i], .7) + COV.get(els[j], .7)):
                s.add((i, j))
    return s


def _integrity(els, R_ref, R_now):
    b_ref, b_now = _bonds(els, R_ref), _bonds(els, R_now)
    broke  = sorted(b for b in (b_ref - b_now) if b not in REACTING)
    formed = sorted(b for b in (b_now - b_ref) if b not in REACTING)
    return ("VALID" if not broke and not formed else "DISTORTED"), broke, formed


def evaluate(name, endpoint, workdir=".", inputs_dir="../inputs"):
    import openmm
    from ash import Fragment, ORCATheory, OpenMMTheory, QMMMTheory, Optimizer

    os.chdir(workdir)
    geomfile = {"reactant": "reactant.xyz", "product": "product.xyz", "ts": "ts.xyz"}[endpoint]
    sub_lines = _read_body(os.path.join(inputs_dir, geomfile), 2, 26)
    els_ref = [l.split()[0] for l in sub_lines]
    R_ref = [[float(v) for v in l.split()[1:4]] for l in sub_lines]

    cd = open("%s_coords.xyz" % name).read().splitlines()
    n = int(cd[0].split()[0])
    site_body = cd[2 + 24:2 + n]
    coords = "\n".join([l.rstrip() for l in sub_lines] + [l.rstrip() for l in site_body])
    mmcharges = [float(x) for x in open("%s_charges.txt" % name).read().split(",")]

    frag = Fragment(coordsstring=coords, charge=-2, mult=1)
    els = [a.split()[0] for a in coords.strip().split("\n")]
    assert len(mmcharges) == frag.numatoms, \
        "charges (%d) != fragment atoms (%d)" % (len(mmcharges), frag.numatoms)

    orca = ORCATheory(orcadir=ORCADIR, orcasimpleinput=SIMPLE, orcablocks=BLOCKS,
                      numcores=NCORES, autostart=False)
    mm = OpenMMTheory(fragment=frag, dummysystem=True, platform="CPU", numcores=1,
                      autoconstraints=None, rigidwater=False)
    lj = openmm.CustomBondForce("4*eps*((sig/r)^12 - (sig/r)^6)")
    lj.addPerBondParameter("sig"); lj.addPerBondParameter("eps"); lj.setUsesPeriodicBoundaryConditions(False)
    ang = openmm.unit.angstrom; kcal = openmm.unit.kilocalorie_per_mole
    for i in range(24):
        si, ei = LJ[els[i]]
        for j in range(24, frag.numatoms):
            sj, ej = LJ[els[j]]
            sig = ((si + sj) / 2 * ang).value_in_unit(openmm.unit.nanometer)
            eps = (((ei * ej) ** 0.5) * kcal).value_in_unit(openmm.unit.kilojoule_per_mole)
            lj.addBond(i, j, [sig, eps])
    lj.setForceGroup(11); mm.system.addForce(lj)

    qmmm = QMMMTheory(qm_theory=orca, mm_theory=mm, fragment=frag, qmatoms=list(range(24)),
                      charges=mmcharges, embedding="elstat", qm_charge=-2, qm_mult=1, numcores=NCORES)
    FROZEN = list(range(24, frag.numatoms))
    Optimizer(theory=qmmm, fragment=frag, coordsystem="hdlc", frozenatoms=FROZEN,
              maxiter=400, conv_criteria=CONV, charge=-2, mult=1)

    R_now = [list(frag.coords[i]) for i in range(24)]
    frag.write_xyzfile("%s_%s_relaxed.xyz" % (name, endpoint))
    verdict, broke, formed = _integrity(els_ref, R_ref, R_now)
    E = float(frag.energy)
    open("%s_%s_result.txt" % (name, endpoint), "w").write(
        "energy_Eh=%.8f\nintegrity=%s\nbroke=%s\nformed=%s\n" % (E, verdict, broke, formed))
    print(">> evaluate %s %s: E=%.6f Eh  integrity=%s  broke=%s formed=%s"
          % (name, endpoint, E, verdict, broke, formed))
    return {"energy": E, "integrity": verdict, "broke": broke, "formed": formed,
            "geometry": list(zip(els_ref, R_now))}


def barrier(*a, **k):
    raise NotImplementedError("NEB barrier belongs to Layer 4.")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: python3 oracle.py <design_name> <reactant|ts|product> [inputs_dir]")
    inp = sys.argv[3] if len(sys.argv) > 3 else "../inputs"
    evaluate(sys.argv[1], sys.argv[2], workdir=".", inputs_dir=inp)
