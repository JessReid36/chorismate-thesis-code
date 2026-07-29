# 8A barrier: is the best intact discrete design actually catalytic? (tests if window metric hides catalysis)

## Protocol (B2: both endpoints under the SAME field, per-design reference)
barrier dE! = E(TS relaxed under 8A field) - E(reactant relaxed under 8A field)
compare to BARE barrier +17.47. If dE! < 17.47 -> catalytic despite being out-of-window (window metric
was hiding it). If dE! >= 17.47 -> confirmed non-catalytic, window verdict stands.

## Steps
1. reactant-under-8A energy: ALREADY DONE (shell8p0 ladder job). Get it:
     bash extract_react_energy.sh   -> E_react (last FINAL SINGLE POINT ENERGY, Eh)
2. build TS-under-8A design + relax it:
     # in a dir with shell8p0_coords.xyz, shell8p0_charges.txt, ts.xyz:
     python3 make_ts_design.py        -> shell8p0ts_coords.xyz + _charges.txt
     python3 relax_barrier.py shell8p0ts     (dry-run)
     qsub run_barrier8_ts.pbs                 -> shell8p0ts_energy.txt (E_TS)
3. barrier = (E_TS - E_react) * 627.509 kcal/mol. Compare to 17.47.

## Then repeat for the fine-sweep winners (7.0-9.0 in 0.1A) once those geometry relaxations are done -
## same make_ts_design + relax_barrier, per shell. Only the shells that stay INTACT get a barrier.
