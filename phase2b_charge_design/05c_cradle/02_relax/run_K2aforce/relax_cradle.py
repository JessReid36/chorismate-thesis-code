import sys, openmm
from ash import Fragment, ORCATheory, OpenMMTheory, QMMMTheory, Optimizer

# 05c cradle Step 2: relax (fixed K-charge design + NEUTRAL 6-bead steric cradle) and measure C1-C6.
# Tests whether a neutral steric pocket pulls the field-opened substrate into the near-attack window
# [2.6-3.5 A] where charges alone could not (K2_aforce 4.06, K6_apot 4.24). Beads are q=0 LJ walls.
#
# Layout: 0-23 QM substrate (-2); then charge surrogates (guan/formate, MM partial charges + element LJ);
# then cradle beads (element 'Ar', q=0, PER-BEAD sigma from *_beadsigma.txt, eps=0.086). All MM frozen.
#
# Usage: relax_cradle.py <name>        dry run ;  relax_cradle.py <name> run   real
#   reads <name>_cradle_coords.xyz, _cradle_charges.txt, _cradle_beadsigma.txt

name=sys.argv[1]; RUN=(len(sys.argv)>2 and sys.argv[2]=="run")
ORCADIR="/home/apps2/ORCA/6.0.1"
coords="\n".join(open("%s_cradle_coords.xyz"%name).read().splitlines()[2:])
mmcharges=[float(x) for x in open("%s_cradle_charges.txt"%name).read().split(",")]
beadsigma=[float(x) for x in open("%s_cradle_beadsigma.txt"%name).read().split(",")]
frag=Fragment(coordsstring=coords, charge=-2, mult=1)
els=[a.split()[0] for a in coords.strip().split("\n")]
nbead=els.count("Ar")
nmm=frag.numatoms-24
print(">> %s: %d atoms (24 QM + %d MM: %d charge-surr + %d beads), MM net %+.2f, QM -2"%(
    name,frag.numatoms,nmm,nmm-nbead,nbead,sum(mmcharges)))
assert len(mmcharges)==frag.numatoms
assert len(beadsigma)==nbead, "beadsigma count %d != nbead %d"%(len(beadsigma),nbead)

qmatoms=list(range(24))
# element LJ (Angstrom sigma, kcal eps) for the CHARGE-surrogate atoms; beads use per-bead sigma
LJ={"C":(3.40,0.086),"N":(3.25,0.170),"O":(2.96,0.210),"H":(1.06,0.016)}
BEAD_EPS=0.086

orca=ORCATheory(orcadir=ORCADIR,
    orcasimpleinput="! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF",
    orcablocks="%cpcm epsilon 4.0 end\n%scf MaxIter 300 Guess PModel end",
    numcores=8, autostart=False)
mm=OpenMMTheory(fragment=frag, dummysystem=True, platform="CPU", numcores=1,
                autoconstraints=None, rigidwater=False)

lj=openmm.CustomBondForce("4*eps*((sig/r)^12 - (sig/r)^6)")
lj.addPerBondParameter("sig"); lj.addPerBondParameter("eps"); lj.setUsesPeriodicBoundaryConditions(False)
ang=openmm.unit.angstrom; kcal=openmm.unit.kilocalorie_per_mole; npair=0
bead_idx=0
for j in range(24,frag.numatoms):
    ej_el=els[j]
    if ej_el=="Ar":
        sj=beadsigma[bead_idx]; ej=BEAD_EPS; bead_idx+=1
    else:
        sj,ej=LJ[ej_el]
    for i in range(24):
        si,ei=LJ[els[i]]
        sig=((si+sj)/2*ang).value_in_unit(openmm.unit.nanometer)
        eps=(((ei*ej)**0.5)*kcal).value_in_unit(openmm.unit.kilojoule_per_mole)
        lj.addBond(i,j,[sig,eps]); npair+=1
lj.setForceGroup(11); mm.system.addForce(lj)
print(">> QM-MM LJ pairs: %d (incl %d bead pairs)"%(npair,nbead*24))

qmmm=QMMMTheory(qm_theory=orca, mm_theory=mm, fragment=frag, qmatoms=qmatoms,
                charges=mmcharges, embedding="elstat", qm_charge=-2, qm_mult=1, numcores=8)
FROZEN=list(range(24,frag.numatoms))
CONV={"convergence_energy":1e-5,"convergence_grms":3.0e-3,"convergence_gmax":6.0e-3,
      "convergence_drms":1.0e-2,"convergence_dmax":1.5e-2}
if not RUN:
    print("\n=== DRY RUN OK [cradle] — %s: %d LJ pairs (%d beads), QM -2, MM %+.0f, %d frozen. ==="%(
        name,npair,nbead,sum(mmcharges),len(FROZEN)))
    sys.exit(0)
Optimizer(theory=qmmm, fragment=frag, coordsystem="hdlc", frozenatoms=FROZEN,
          maxiter=400, conv_criteria=CONV, charge=-2, mult=1)
frag.write_xyzfile("%s_cradle_relaxed.xyz"%name)
print(">> DONE -> %s_cradle_relaxed.xyz"%name)
