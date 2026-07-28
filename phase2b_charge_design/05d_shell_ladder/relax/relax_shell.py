import sys, openmm
from ash import Fragment, ORCATheory, OpenMMTheory, QMMMTheory, Optimizer
name=sys.argv[1]; RUN=(len(sys.argv)>2 and sys.argv[2]=="run")
ORCADIR="/home/apps2/ORCA/6.0.1"
coords="\n".join(open("%s_coords.xyz"%name).read().splitlines()[2:])
mmcharges=[float(x) for x in open("%s_charges.txt"%name).read().split(",")]
frag=Fragment(coordsstring=coords, charge=-2, mult=1)
els=[a.split()[0] for a in coords.strip().split("\n")]
print(">> %s: %d atoms (24 QM + %d MM), MM net %+.2f"%(name,frag.numatoms,frag.numatoms-24,sum(mmcharges)))
assert len(mmcharges)==frag.numatoms
LJ={"C":(3.40,0.086),"N":(3.25,0.170),"O":(2.96,0.210),"H":(1.06,0.016)}
orca=ORCATheory(orcadir=ORCADIR,orcasimpleinput="! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF",
    orcablocks="%cpcm epsilon 4.0 end\n%scf MaxIter 300 Guess PModel end",numcores=8,autostart=False)
mm=OpenMMTheory(fragment=frag,dummysystem=True,platform="CPU",numcores=1,autoconstraints=None,rigidwater=False)
lj=openmm.CustomBondForce("4*eps*((sig/r)^12 - (sig/r)^6)")
lj.addPerBondParameter("sig");lj.addPerBondParameter("eps");lj.setUsesPeriodicBoundaryConditions(False)
ang=openmm.unit.angstrom;kcal=openmm.unit.kilocalorie_per_mole;npair=0
for i in range(24):
    si,ei=LJ[els[i]]
    for j in range(24,frag.numatoms):
        sj,ej=LJ[els[j]]
        sig=((si+sj)/2*ang).value_in_unit(openmm.unit.nanometer)
        eps=(((ei*ej)**0.5)*kcal).value_in_unit(openmm.unit.kilojoule_per_mole)
        lj.addBond(i,j,[sig,eps]);npair+=1
lj.setForceGroup(11);mm.system.addForce(lj)
print(">> LJ pairs: %d"%npair)
qmmm=QMMMTheory(qm_theory=orca,mm_theory=mm,fragment=frag,qmatoms=list(range(24)),charges=mmcharges,
    embedding="elstat",qm_charge=-2,qm_mult=1,numcores=8)
FROZEN=list(range(24,frag.numatoms))
CONV={"convergence_energy":1e-5,"convergence_grms":3.0e-3,"convergence_gmax":6.0e-3,"convergence_drms":1.0e-2,"convergence_dmax":1.5e-2}
if not RUN:
    print("\n=== DRY RUN OK %s: %d pairs, MM %+.0f, %d frozen ==="%(name,npair,sum(mmcharges),len(FROZEN))); sys.exit(0)
Optimizer(theory=qmmm,fragment=frag,coordsystem="hdlc",frozenatoms=FROZEN,maxiter=400,conv_criteria=CONV,charge=-2,mult=1)
frag.write_xyzfile("%s_relaxed.xyz"%name); print(">> DONE -> %s_relaxed.xyz"%name)
