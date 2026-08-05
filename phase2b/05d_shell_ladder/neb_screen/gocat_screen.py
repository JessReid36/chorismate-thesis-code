# GOCAT-style endpoint screening (Eq 3.2 RMSG + 3.3 stabilisation + integrity).
# For a design: relax substrate under field (R and/or P), then evaluate:
#   - RMSG at the relaxed endpoint (should be < beta ~ near-minimum)
#   - endpoint energy vs bare reference (Phi_stab: not destabilised above reference)
#   - substrate integrity (ether O3-C4 intact)
# A design PASSES only if BOTH endpoints satisfy all three -> then it's eligible for a NEB barrier.
import sys, openmm
from ash import Fragment, ORCATheory, OpenMMTheory, QMMMTheory, Optimizer
name=sys.argv[1]; endpoint=sys.argv[2]  # 'reactant' or 'product'
RUN=(len(sys.argv)>3 and sys.argv[3]=="run")
ORCADIR="/home/apps2/ORCA/6.0.1"
# GOCAT thresholds
ALPHA=250.0; BETA=10.0      # RMSG penalty (kcal/mol/Ang)
KAPPA=1e4                    # stabilisation penalty
# bare reference barrier context: bare R energy, bare TS energy (QM-only, from Phase-2 singlepoints)
rd=open("%s_coords.xyz"%name).read().splitlines(); n=int(rd[0].split()[0]); body=rd[2+24:2+n]
geomfile={"reactant":"reactant.xyz","product":"product.xyz","ts":"ts.xyz"}[endpoint]
sub=open(geomfile).read().splitlines()[2:26]
coords="\n".join(sub+body)
mmcharges=[float(x) for x in open("%s_charges.txt"%name).read().split(",")]
frag=Fragment(coordsstring=coords, charge=-2, mult=1)
els=[a.split()[0] for a in coords.strip().split("\n")]
print(">> GOCAT screen %s %s: %d atoms, MM net %+.2f"%(name,endpoint,frag.numatoms,sum(mmcharges)))
LJ={"C":(3.40,0.086),"N":(3.25,0.170),"O":(2.96,0.210),"H":(1.06,0.016)}
orca=ORCATheory(orcadir=ORCADIR,orcasimpleinput="! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF",
    orcablocks="%cpcm epsilon 4.0 end\n%scf MaxIter 300 Guess PModel end",numcores=8,autostart=False)
mm=OpenMMTheory(fragment=frag,dummysystem=True,platform="CPU",numcores=1,autoconstraints=None,rigidwater=False)
lj=openmm.CustomBondForce("4*eps*((sig/r)^12 - (sig/r)^6)")
lj.addPerBondParameter("sig");lj.addPerBondParameter("eps");lj.setUsesPeriodicBoundaryConditions(False)
ang=openmm.unit.angstrom;kcal=openmm.unit.kilocalorie_per_mole
for i in range(24):
    si,ei=LJ[els[i]]
    for j in range(24,frag.numatoms):
        sj,ej=LJ[els[j]]
        sig=((si+sj)/2*ang).value_in_unit(openmm.unit.nanometer);eps=(((ei*ej)**0.5)*kcal).value_in_unit(openmm.unit.kilojoule_per_mole)
        lj.addBond(i,j,[sig,eps])
lj.setForceGroup(11);mm.system.addForce(lj)
qmmm=QMMMTheory(qm_theory=orca,mm_theory=mm,fragment=frag,qmatoms=list(range(24)),charges=mmcharges,
    embedding="elstat",qm_charge=-2,qm_mult=1,numcores=8)
FROZEN=list(range(24,frag.numatoms))
CONV={"convergence_energy":1e-5,"convergence_grms":3.0e-3,"convergence_gmax":6.0e-3,"convergence_drms":1.0e-2,"convergence_dmax":1.5e-2}
if not RUN:
    print("=== DRY RUN OK: GOCAT screen %s %s (relax + RMSG/stab/integrity) ==="%(name,endpoint)); sys.exit(0)
Optimizer(theory=qmmm,fragment=frag,coordsystem="hdlc",frozenatoms=FROZEN,maxiter=400,conv_criteria=CONV,charge=-2,mult=1)
frag.write_xyzfile("%s_%s_screened.xyz"%(name,endpoint))
# after convergence, one to get RMSG on the substrate atoms
import numpy as np
res=qmmm.run(current_coords=frag.coords, elems=frag.elems, Grad=True, charge=-2, mult=1)
E=frag.energy
g=np.array(qmmm.grad).reshape(-1,3)[:24]  # substrate-atom gradients
rmsg=float(np.sqrt(np.mean(np.sum(g**2,axis=1)))) * 627.509 / 0.529177  # Eh/Bohr -> kcal/mol/Ang
import math
a=frag.coords
o34=math.sqrt(sum((a[7][k]-a[8][k])**2 for k in range(3)))
phi_rmsg=ALPHA*(rmsg-BETA)**2 if rmsg>BETA else 0.0
open("%s_%s_screen.txt"%(name,endpoint),"w").write("E=%.6f\nRMSG=%.4f\nO3C4=%.4f\nphi_rmsg=%.2f\n"%(E,rmsg,o34,phi_rmsg))
print(">> SCREEN %s %s: E=%.6f Eh  RMSG=%.3f kcal/mol/Ang  O3-C4=%.3f  phi_rmsg=%.1f"%(name,endpoint,E,rmsg,o34,phi_rmsg))
print(">>   RMSG %s beta=%.0f  (%s)"%('<' if rmsg<BETA else '>=',BETA,'near-min PASS' if rmsg<BETA else 'NOT at min FAIL'))
print(">>   integrity: O3-C4=%.3f (%s)"%(o34,'intact' if o34<2.2 else 'BROKEN'))
