import sys, openmm
from ash import Fragment, ORCATheory, OpenMMTheory, QMMMTheory, Optimizer

RUN = (len(sys.argv) > 1 and sys.argv[1] == "run")
ORCADIR = "/home/apps2/ORCA/6.0.1"
SIGMA_A = 3.8
EPS     = 1.0

sub = """C  60.11193800 34.38396900 57.21013200
H  60.60719900 33.42993900 57.38460900
H  60.72469300 35.27181800 57.05919300
C  58.78057700 34.50353200 57.15195800
C  58.09218200 35.83426100 56.89335400
O  58.81246700 36.81821800 56.57878900
O  56.83954900 35.84891700 57.00961700
O  57.93461500 33.42640300 57.34599800
C  57.67693600 32.53773200 56.20080600
H  57.01901700 31.76376900 56.62407200
C  58.98648700 31.96184000 55.79317200
H  59.46393000 31.24298300 56.46049600
C  59.68310500 32.48267100 54.76815100
C  59.06432100 33.48370900 53.89317900
H  59.70236100 33.90454400 53.11096900
C  57.78125300 33.86565900 54.03153300
H  57.33218100 34.60235400 53.35662900
C  56.88979900 33.27109300 55.08618200
H  56.34719000 34.09064800 55.58588000
O  55.96286300 32.39770700 54.43466800
H  55.31084500 32.10952700 55.09277000
C  61.12320000 32.06212600 54.53420800
O  61.78496200 31.70232800 55.54338100
O  61.54107100 32.10651300 53.35273800"""
charge_sites = [(-1.0, 63.004400, 29.116300, 52.225600),
                ( 1.0, 57.499400, 32.313000, 60.999500)]

lines = sub.strip().split("\n")
for q,x,y,z in charge_sites:
    lines.append("He %.6f %.6f %.6f" % (x,y,z))
frag = Fragment(coordsstring="\n".join(lines), charge=-2, mult=1)
print(">> Fragment: %d atoms, QM charge -2" % frag.numatoms)

qmatoms = list(range(24))
mmcharges = [0.0]*24 + [q for (q,_,_,_) in charge_sites]

orca = ORCATheory(orcadir=ORCADIR,
                  orcasimpleinput="! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF",
                  orcablocks="%cpcm epsilon 4.0 end\n%scf MaxIter 300 Guess PModel end",
                  numcores=8, autostart=False)
print(">> ORCATheory built")

mm = OpenMMTheory(fragment=frag, dummysystem=True, platform="CPU", numcores=1,
                  autoconstraints=None, rigidwater=False)
print(">> OpenMMTheory built")

wall = openmm.CustomBondForce("eps*(sigma/r)^12")
wall.addPerBondParameter("sigma"); wall.addPerBondParameter("eps")
wall.setUsesPeriodicBoundaryConditions(False)
signm = (SIGMA_A*openmm.unit.angstrom).value_in_unit(openmm.unit.nanometer)
epskj = (EPS*openmm.unit.kilocalorie_per_mole).value_in_unit(openmm.unit.kilojoule_per_mole)
npair=0
for a in range(24):
    for c in (24,25):
        wall.addBond(a, c, [signm, epskj]); npair+=1
wall.setForceGroup(11)
mm.system.addForce(wall)
print(">> smooth wall added: %d pairs, sigma=%.1f A r^-12 no-step" % (npair, SIGMA_A))

qmmm = QMMMTheory(qm_theory=orca, mm_theory=mm, fragment=frag,
                  qmatoms=qmatoms, charges=mmcharges, embedding="elstat",
                  qm_charge=-2, qm_mult=1, numcores=8)
print(">> QMMMTheory built | QM stays -2")

if not RUN:
    print("\n=== DRY RUN OK (smooth-wall) — %d pairs sigma %.1f A. ===" % (npair, SIGMA_A))
    sys.exit(0)

Optimizer(theory=qmmm, fragment=frag, coordsystem="hdlc",
          maxiter=200, charge=-2, mult=1)
frag.write_xyz("k2_wall3_relaxed.xyz")
print(">> DONE -> k2_wall3_relaxed.xyz")
