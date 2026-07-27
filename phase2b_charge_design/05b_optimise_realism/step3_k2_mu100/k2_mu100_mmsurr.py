import sys, openmm
from ash import Fragment, ORCATheory, OpenMMTheory, QMMMTheory, Optimizer

# Step 3 (PLAN.md) — RELAXED VALIDATION of the K2 mu=100 PREORGANISED design (the arbiter).
# Design (05b Step 2, certified gap 0.000): TWO +1 at sites 105/121 (net surrogate +2), D=+0.0006 (~0
# distorting drive) vs the mu=0 distorting control (sites 150/249, D=+0.105) which relaxed to C1-C6=5.0 A.
# Question: does pricing in preorganisation (a_pot penalty) make the relaxed C1-C6 stay COMPACT (2.6-3.5 A)?
# If yes -> the fix works, a_pot was the right mode. If it still opens -> revisit with a_force.
#
# Method = the validated MM-molecular-surrogate production method (k2_mmsurr_spike.py), adapted:
#   44 atoms: 0-23 QM substrate (dianion -2); 24-33 MM guanidinium(+1) @site105; 34-43 MM guanidinium(+1)
#   @site121. NO formate (design is two +1). QM stays -2 (MM adds field, not electrons). Each guanidinium:
#   partial charges (C+0.64, per-N -0.80, per-H +0.46) AND explicit per-pair LJ vs every QM atom.
# Steric realizability CONFIRMED: guan-guan min 2.85 A, guan-substrate min 2.20-2.38 A (no clash).

RUN = (len(sys.argv) > 1 and sys.argv[1] == "run")
ORCADIR = "/home/apps2/ORCA/6.0.1"

coords = """C     60.111938    34.383969    57.210132
H     60.607199    33.429939    57.384609
H     60.724693    35.271818    57.059193
C     58.780577    34.503532    57.151958
C     58.092182    35.834261    56.893354
O     58.812467    36.818218    56.578789
O     56.839549    35.848917    57.009617
O     57.934615    33.426403    57.345998
C     57.676936    32.537732    56.200806
H     57.019017    31.763769    56.624072
C     58.986487    31.961840    55.793172
H     59.463930    31.242983    56.460496
C     59.683105    32.482671    54.768151
C     59.064321    33.483709    53.893179
H     59.702361    33.904544    53.110969
C     57.781253    33.865659    54.031533
H     57.332181    34.602354    53.356629
C     56.889799    33.271093    55.086182
H     56.347190    34.090648    55.585880
O     55.962863    32.397707    54.434668
H     55.310845    32.109527    55.092770
C     61.123200    32.062126    54.534208
O     61.784962    31.702328    55.543381
O     61.541071    32.106513    53.352738
C     54.837300    30.677200    58.906400
N     54.987877    31.062420    57.642346
H     54.677360    31.724922    58.338600
H     55.184046    30.107382    57.906013
N     54.428401    31.549604    59.823249
H     54.624570    30.594566    60.086916
H     54.542749    31.842140    58.863328
N     55.095621    29.419576    59.253605
H     55.209970    29.712111    58.293684
H     54.785104    30.082078    59.949859
C     57.971400    29.898500    59.181500
N     58.121977    30.283720    57.917446
H     57.811460    30.946222    58.613700
H     58.318146    29.328682    58.181113
N     57.562501    30.770904    60.098349
H     57.758670    29.815866    60.362016
H     57.676849    31.063440    59.138428
N     58.229721    28.640876    59.528705
H     58.344070    28.933411    58.568784
H     57.919204    29.303378    60.224959"""

frag = Fragment(coordsstring=coords, charge=-2, mult=1)
print(">> Fragment: %d atoms (24 QM + 20 MM = 2 guanidinium)" % frag.numatoms)

qmatoms = list(range(24))
guan_q = [0.64] + [-0.80,0.46,0.46]*3          # net +1 each
mmcharges = [0.0]*24 + guan_q + guan_q          # two guanidiniums -> MM net +2
print(">> MM charge sum = %+.2f (net +2 = two +1 guanidiniums); QM subsystem stays -2" % sum(mmcharges))

LJ = {"C":(3.40,0.086),"N":(3.25,0.170),"O":(2.96,0.210),"H":(1.06,0.016)}
els = [a.split()[0] for a in coords.strip().split("\n")]

orca = ORCATheory(orcadir=ORCADIR,
                  orcasimpleinput="! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF",
                  orcablocks="%cpcm epsilon 4.0 end\n%scf MaxIter 300 Guess PModel end",
                  numcores=8, autostart=False)
print(">> ORCATheory built")

mm = OpenMMTheory(fragment=frag, dummysystem=True, platform="CPU", numcores=1,
                  autoconstraints=None, rigidwater=False)
print(">> OpenMMTheory built")

# explicit QM-substrate <-> MM-surrogate Lennard-Jones (LB combining): wall + directional H-bond well
lj = openmm.CustomBondForce("4*eps*((sig/r)^12 - (sig/r)^6)")
lj.addPerBondParameter("sig"); lj.addPerBondParameter("eps")
lj.setUsesPeriodicBoundaryConditions(False)
ang = openmm.unit.angstrom; kcal = openmm.unit.kilocalorie_per_mole
npair=0
for i in range(24):
    si,ei = LJ[els[i]]
    for j in range(24,44):
        sj,ej = LJ[els[j]]
        sig = ((si+sj)/2*ang).value_in_unit(openmm.unit.nanometer)
        eps = (((ei*ej)**0.5)*kcal).value_in_unit(openmm.unit.kilojoule_per_mole)
        lj.addBond(i, j, [sig, eps]); npair+=1
lj.setForceGroup(11)
mm.system.addForce(lj)
print(">> explicit QM-MM LJ added: %d pairs" % npair)

qmmm = QMMMTheory(qm_theory=orca, mm_theory=mm, fragment=frag,
                  qmatoms=qmatoms, charges=mmcharges, embedding="elstat",
                  qm_charge=-2, qm_mult=1, numcores=8)
print(">> QMMMTheory built | QM -2 | MM two-guanidinium surrogate (charges+LJ)")

FROZEN = list(range(24,44))  # both guanidiniums fixed at certified positions

# LOOSE conv (floppy field-distorted surface; validated for the surrogate method)
CONV = {"convergence_energy":1e-5,"convergence_grms":3.0e-3,"convergence_gmax":6.0e-3,
        "convergence_drms":1.0e-2,"convergence_dmax":1.5e-2}

if not RUN:
    print("\n=== DRY RUN OK — %d LJ pairs, QM -2, MM +2, %d frozen. ===" % (npair,len(FROZEN)))
    sys.exit(0)

Optimizer(theory=qmmm, fragment=frag, coordsystem="hdlc", frozenatoms=FROZEN,
          maxiter=300, conv_criteria=CONV, charge=-2, mult=1)
frag.write_xyzfile("k2_mu100_relaxed.xyz")
print(">> DONE -> k2_mu100_relaxed.xyz")
