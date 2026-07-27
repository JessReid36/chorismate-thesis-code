import sys, openmm
from ash import Fragment, ORCATheory, OpenMMTheory, QMMMTheory, Optimizer

# Step 3 (a_force fallback) — RELAXED VALIDATION of the a_force-preorganised K2 design.
# After a_pot mu=100 FAILED (C1-C6 relaxed to 5.28 A, ~same as mu=0 control), test the force-imbalance
# coordinate. a_force mu=100 keeps a DIPOLE: guanidinium(+1)@249 (unchanged from orig K2) + formate(-1)
# relocated 150->185. D_force=0.006 (certified gap 0.000). Net-neutral MM (dipole), QM stays -2.
# Steric realizability CONFIRMED (guan-formate 9.33, guan/formate-substrate 2.51/2.59 A).
#
# DECISION: if this relaxes C1-C6 into [2.6,3.5] -> a_force is the right preorg mode, FIX WORKS.
#           if it ALSO stays open -> no K=2 monopole pair can preorganise this dianion -> higher-K/cradle.
#
# Method = validated MM-molecular-surrogate production (k2_mmsurr_spike.py), guanidinium+formate.

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
C     57.499400    32.313000    60.999500
N     57.649977    32.698220    59.735446
H     57.339460    33.360722    60.431700
H     57.846146    31.743182    59.999113
N     57.090501    33.185404    61.916349
H     57.286670    32.230366    62.180016
H     57.204849    33.477940    60.956428
N     57.757721    31.055376    61.346705
H     57.872070    31.347911    60.386784
H     57.447204    31.717878    62.042959
C     64.308600    30.099700    52.535500
O     65.492425    30.161476    52.138963
O     63.597316    29.072317    52.568059
H     63.850621    31.035550    52.888261"""

frag = Fragment(coordsstring=coords, charge=-2, mult=1)
print(">> Fragment: %d atoms (24 QM + 10 guan + 4 formate)" % frag.numatoms)

qmatoms = list(range(24))
guan_q = [0.64] + [-0.80,0.46,0.46]*3            # +1
form_q = [0.45, -0.80, -0.80, 0.15]              # -1 (C,O,O,H)
mmcharges = [0.0]*24 + guan_q + form_q           # MM net 0 (dipole)
print(">> MM charge sum = %+.2f (net 0 dipole); QM stays -2" % sum(mmcharges))

LJ = {"C":(3.40,0.086),"N":(3.25,0.170),"O":(2.96,0.210),"H":(1.06,0.016)}
els = [a.split()[0] for a in coords.strip().split("\n")]

orca = ORCATheory(orcadir=ORCADIR,
                  orcasimpleinput="! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF",
                  orcablocks="%cpcm epsilon 4.0 end\n%scf MaxIter 300 Guess PModel end",
                  numcores=8, autostart=False)
mm = OpenMMTheory(fragment=frag, dummysystem=True, platform="CPU", numcores=1,
                  autoconstraints=None, rigidwater=False)

lj = openmm.CustomBondForce("4*eps*((sig/r)^12 - (sig/r)^6)")
lj.addPerBondParameter("sig"); lj.addPerBondParameter("eps")
lj.setUsesPeriodicBoundaryConditions(False)
ang = openmm.unit.angstrom; kcal = openmm.unit.kilocalorie_per_mole
npair=0
for i in range(24):
    si,ei = LJ[els[i]]
    for j in range(24,38):
        sj,ej = LJ[els[j]]
        sig = ((si+sj)/2*ang).value_in_unit(openmm.unit.nanometer)
        eps = (((ei*ej)**0.5)*kcal).value_in_unit(openmm.unit.kilojoule_per_mole)
        lj.addBond(i, j, [sig, eps]); npair+=1
lj.setForceGroup(11); mm.system.addForce(lj)
print(">> QM-MM LJ pairs: %d" % npair)

qmmm = QMMMTheory(qm_theory=orca, mm_theory=mm, fragment=frag, qmatoms=qmatoms,
                  charges=mmcharges, embedding="elstat", qm_charge=-2, qm_mult=1, numcores=8)
FROZEN = list(range(24,38))
CONV = {"convergence_energy":1e-5,"convergence_grms":3.0e-3,"convergence_gmax":6.0e-3,
        "convergence_drms":1.0e-2,"convergence_dmax":1.5e-2}

if not RUN:
    print("\n=== DRY RUN OK — %d LJ pairs, QM -2, MM 0 (dipole), %d frozen. ===" % (npair,len(FROZEN)))
    sys.exit(0)

Optimizer(theory=qmmm, fragment=frag, coordsystem="hdlc", frozenatoms=FROZEN,
          maxiter=300, conv_criteria=CONV, charge=-2, mult=1)
frag.write_xyzfile("k2_mu100_aforce_relaxed.xyz")
print(">> DONE -> k2_mu100_aforce_relaxed.xyz")
