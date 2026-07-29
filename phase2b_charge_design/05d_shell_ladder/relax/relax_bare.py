import sys
from ash import Fragment, ORCATheory, Optimizer
ORCADIR="/home/apps2/ORCA/6.0.1"
RUN=(len(sys.argv)>1 and sys.argv[1]=="run")
# bare 24-atom substrate reactant geometry
coords="\n".join(open("reactant.xyz").read().splitlines()[2:26])
frag=Fragment(coordsstring=coords, charge=-2, mult=1)
print(">> BARE control: %d atoms (should be 24), charge -2"%frag.numatoms)
assert frag.numatoms==24
orca=ORCATheory(orcadir=ORCADIR,orcasimpleinput="! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF",
    orcablocks="%cpcm epsilon 4.0 end\n%scf MaxIter 300 Guess PModel end",numcores=8,autostart=False)
CONV={"convergence_energy":1e-5,"convergence_grms":3.0e-3,"convergence_gmax":6.0e-3,"convergence_drms":1.0e-2,"convergence_dmax":1.5e-2}
if not RUN:
    print("=== DRY RUN OK: bare 24-atom substrate, pure QM (no bodies, no field) ==="); sys.exit(0)
# NO frozen atoms (whole substrate free), pure QM optimisation
Optimizer(theory=orca,fragment=frag,coordsystem="hdlc",maxiter=400,conv_criteria=CONV,charge=-2,mult=1)
frag.write_xyzfile("bare_control_relaxed.xyz"); print(">> DONE -> bare_control_relaxed.xyz")
