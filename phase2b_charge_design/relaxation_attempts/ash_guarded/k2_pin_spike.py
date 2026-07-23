import sys, os
from ash import Fragment, ORCATheory, OpenMMTheory, QMMMTheory, Optimizer

RUN = (len(sys.argv) > 1 and sys.argv[1] == "run")
ORCADIR = "/home/apps2/ORCA/6.0.1"

# 24-atom substrate (reactant), 0-based indices; O3=7, C4=8, C1=0, C6=12
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

# K2 charge sites appended AFTER the 24 substrate atoms:
#   idx 24 = -1 at 8.40 A ; idx 25 = +1 at 3.84 A (the Arg90 valley, the implosion driver)
charge_sites = [(-1.0, 63.004400, 29.116300, 52.225600),
                ( 1.0, 57.499400, 32.313000, 60.999500)]

# build combined coordinate string: substrate as real atoms, charge sites as He placeholders
# (element is irrelevant for MM point charges; QM region is only atoms 0-23)
lines = sub.strip().split("\n")
for q,x,y,z in charge_sites:
    lines.append("He %.6f %.6f %.6f" % (x,y,z))
allcoords = "\n".join(lines)

frag = Fragment(coordsstring=allcoords, charge=-2, mult=1)
print(">> Fragment built: %d atoms total (24 QM + %d MM charge sites)" % (frag.numatoms, len(charge_sites)))

qmatoms = list(range(24))
mmcharges = [0.0]*24 + [q for (q,_,_,_) in charge_sites]   # QM atoms carry 0 in the MM charge list
print(">> qmatoms:", qmatoms[0], "..", qmatoms[-1], "| MM charges on sites:", mmcharges[24:])

orca = ORCATheory(orcadir=ORCADIR,
                  orcasimpleinput="! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF",
                  orcablocks="%cpcm epsilon 4.0 end\n%scf MaxIter 300 end",
                  numcores=8)
print(">> ORCATheory built (level unchanged; CPCM eps=4)")

# minimal MM theory carrying only the charge sites + LJ (sigma/eps ~0 => pure point charges)
mm = OpenMMTheory(fragment=frag, dummysystem=True, platform="CPU", numcores=1, autoconstraints=None, rigidwater=False)
print(">> OpenMMTheory (dummysystem) built")

qmmm = QMMMTheory(qm_theory=orca, mm_theory=mm, fragment=frag,
                  qmatoms=qmatoms, charges=mmcharges, embedding="elstat",
                  qm_charge=-2, qm_mult=1, numcores=8)
print(">> QMMMTheory built | embedding=elstat")

# GUARD: freeze the imploding O3(7)->+1(25) distance at the design value 3.84 A
constraints = {"bond": [[7, 25, 3.84]]}
print(">> constraint:", constraints, "(fixed-value; frame test)")

if not RUN:
    print("\n=== DRY RUN OK — all objects constructed. Re-run with 'run' to optimize. ===")
    sys.exit(0)

Optimizer(theory=qmmm, fragment=frag, coordsystem="hdlc",
          constraints=constraints, constrainvalue=True, maxiter=60,
          charge=-2, mult=1)
frag.write_xyz("k2_guarded_relaxed.xyz")
print(">> DONE. Relaxed geometry written to k2_guarded_relaxed.xyz")
