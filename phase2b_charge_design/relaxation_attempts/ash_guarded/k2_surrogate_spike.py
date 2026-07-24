import sys
from ash import Fragment, ORCATheory, Optimizer

RUN = (len(sys.argv) > 1 and sys.argv[1] == "run")
ORCADIR = "/home/apps2/ORCA/6.0.1"

# 38 atoms: 0-23 substrate (dianion), 24-33 guanidinium(+1), 34-37 formate(-1); net -2, singlet.
# Surrogates are REAL QM molecules at the certified +1/-1 positions -> finite size + Pauli wall.
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
C     63.004400    29.116300    52.225600
O     64.188225    29.178076    51.829063
O     62.293116    28.088917    52.258159
H     62.546421    30.052150    52.578361"""

frag = Fragment(coordsstring=coords, charge=-2, mult=1)
print(">> Fragment: %d atoms (24 substrate + 10 guanidinium + 4 formate), charge -2" % frag.numatoms)

FROZEN = list(range(24, 38))   # freeze all surrogate atoms; relax only the substrate
print(">> frozen surrogate atoms:", FROZEN[0], "..", FROZEN[-1])

orca = ORCATheory(orcadir=ORCADIR,
                  orcasimpleinput="! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM TightSCF",
                  orcablocks="%cpcm epsilon 4.0 end\n%scf MaxIter 300 Guess PModel end",
                  numcores=8, autostart=False)
print(">> ORCATheory built (level unchanged; CPCM eps=4)")

if not RUN:
    print("\n=== DRY RUN OK (surrogate version) — 38-atom QM fragment, 14 frozen. Re-run with 'run'. ===")
    sys.exit(0)

Optimizer(theory=orca, fragment=frag, coordsystem="hdlc",
          frozenatoms=FROZEN, maxiter=80, charge=-2, mult=1)
frag.write_xyz("k2_surrogate_relaxed.xyz")
print(">> DONE. Relaxed geometry -> k2_surrogate_relaxed.xyz")
