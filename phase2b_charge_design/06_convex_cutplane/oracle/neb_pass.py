import sys, openmm
from ash import Fragment, ORCATheory, OpenMMTheory, QMMMTheory, NEBTS
name=sys.argv[1]; RUN=(len(sys.argv)>2 and sys.argv[2]=="run")
ORCADIR="/home/apps2/ORCA/6.0.1"
r_coords="\n".join(open("%s_relaxed.xyz"%name).read().splitlines()[2:])
p_coords="\n".join(open("%s_prod_relaxed.xyz"%name).read().splitlines()[2:])
mmcharges=[float(x) for x in open("%s_charges.txt"%name).read().split(",")]
fragR=Fragment(coordsstring=r_coords, charge=-2, mult=1)
fragP=Fragment(coordsstring=p_coords, charge=-2, mult=1)
els=[a.split()[0] for a in r_coords.strip().split("\n")]
print(">> %s NEB(pass): FIELD-RELAXED endpoints. %d atoms (24 QM + %d MM), MM net %+.2f"%(name,fragR.numatoms,fragR.numatoms-24,sum(mmcharges)))
assert fragR.numatoms==fragP.numatoms==len(mmcharges)
LJ={"C":(3.40,0.086),"N":(3.25,0.170),"O":(2.96,0.210),"H":(1.06,0.016)}
orca=ORCATheory(orcadir=ORCADIR,orcasimpleinput="! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF",
    orcablocks="%cpcm epsilon 4.0 end\n%scf MaxIter 300 Guess PModel end\n%geom coordsys cartesian end",
    numcores=8,autostart=False)
mm=OpenMMTheory(fragment=fragR,dummysystem=True,platform="CPU",numcores=1,autoconstraints=None,rigidwater=False)
lj=openmm.CustomBondForce("4*eps*((sig/r)^12 - (sig/r)^6)")
lj.addPerBondParameter("sig");lj.addPerBondParameter("eps");lj.setUsesPeriodicBoundaryConditions(False)
ang=openmm.unit.angstrom;kcal=openmm.unit.kilocalorie_per_mole
for i in range(24):
    si,ei=LJ[els[i]]
    for j in range(24,fragR.numatoms):
        sj,ej=LJ[els[j]]
        sig=((si+sj)/2*ang).value_in_unit(openmm.unit.nanometer);eps=(((ei*ej)**0.5)*kcal).value_in_unit(openmm.unit.kilojoule_per_mole)
        lj.addBond(i,j,[sig,eps])
lj.setForceGroup(11);mm.system.addForce(lj)
qmmm=QMMMTheory(qm_theory=orca,mm_theory=mm,fragment=fragR,qmatoms=list(range(24)),charges=mmcharges,
    embedding="elstat",qm_charge=-2,qm_mult=1,numcores=8)
if not RUN:
    print("=== DRY RUN OK %s NEB(pass) ==="%name); sys.exit(0)
NEBTS(reactant=fragR, product=fragP, theory=qmmm, images=12, CI=True,
      charge=-2, mult=1, printlevel=2, maxiter=200, free_end=True,
      ActiveRegion=True, actatoms=list(range(24)))
print(">> NEB DONE for %s"%name)
