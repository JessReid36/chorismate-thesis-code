#!/bin/bash
# Stage B: NEB-under-field for the 14 designs that PASSED the GOCAT screen.
# Uses the SCREENED endpoints (converged field-relaxed R + P) as the band endpoints -> NEB will converge.
# Run from 05d_shell_ladder/relax.
set -e
PASS="2p0 7p1 7p3 7p4 7p6 7p8 8p1 8p2 8p5 8p6 8p9 12p0 13p0 15p0"
for tag in $PASS; do
  dtag="shell$tag"
  R=run_screen_${tag}_reactant/${dtag}_reactant_screened.xyz
  P=run_screen_${tag}_product/${dtag}_product_screened.xyz
  if [ ! -f "$R" ] || [ ! -f "$P" ]; then echo "SKIP $tag (missing screened endpoint)"; continue; fi
  # charges
  CH=${dtag}_charges.txt
  [ -f "$CH" ] || CH=run_shell${tag}/${dtag}_charges.txt
  [ -f "$CH" ] || { echo "SKIP $tag (no charges)"; continue; }
  mkdir -p run_nebP_$tag
  cp "$R" run_nebP_$tag/${dtag}_relaxed.xyz          # reactant endpoint (screened)
  cp "$P" run_nebP_$tag/${dtag}_prod_relaxed.xyz     # product endpoint (screened)
  cp "$CH" run_nebP_$tag/${dtag}_charges.txt
  # v2 NEB script (field-relaxed endpoints, free_end, maxiter 200, ActiveRegion freeze)
  { cat <<'PYEOF'
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
PYEOF
  } > run_nebP_$tag/neb_pass.py
  { cat <<PBS
#!/bin/bash
#PBS -N nebP_$tag
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -l walltime=96:00:00
#PBS -M 18660916@sun.ac.za -m ae
#PBS -o nebP_$tag.stdout
#PBS -e nebP_$tag.stderr
cd \$PBS_O_WORKDIR
source /apps/mambaforge/etc/profile.d/conda.sh
conda activate \$HOME/envs/ash
export PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PATH=/apps/openmpi/4.1.1/bin:\$PATH
export LD_LIBRARY_PATH=/home/apps2/ORCA/6.0.1/lib:/apps/openmpi/4.1.1/lib:/apps/mambaforge/envs/medaka/lib:\$LD_LIBRARY_PATH
python3 neb_pass.py $dtag run > nebP_${tag}.log 2>&1
PBS
  } > run_nebP_$tag/run_nebP_$tag.pbs
  ( cd run_nebP_$tag && qsub run_nebP_$tag.pbs )
  echo "submitted nebP_$tag"
done
echo "Stage B (NEB on PASS designs) submitted."
