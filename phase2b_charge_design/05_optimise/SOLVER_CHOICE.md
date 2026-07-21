# Solver choice: HiGHS over CBC for the Tier-1 certified screen

## The requirement
Tier-1's contribution is CERTIFIED optimality -- branch-and-bound gap = 0.000 proves no arrangement
beats the reported one under the objective. This certificate is the differentiator vs GOCAT's genetic
algorithm (a GA has no bound and no gap). A row with an OPEN gap is a row where we have temporarily
lost the one thing that distinguishes the method from a heuristic. So closing the gap is not cosmetic.

## The problem
python-mip's default engine is CBC (COIN-OR). CBC closes the simpler MAX-LOWERING objective at
gap 0.000 for all K, but STALLS on the harder DISTRIBUTED min-max objective (the Mx variable + band
constraints): gaps of 0.04-0.95 remained open at K>=6 even at 120-600 s. The answer was correct --
the gap just would not close.

## The benchmark (solver_benchmark.py, reproducible; result in solver_benchmark_result.txt)
Same distributed min-max MILP, ~330 points, K=7, 120 s budget:
  CBC   : obj 0.7827, gap 0.9519 (NOT certified) at 120 s (hit limit)
  HiGHS : obj 0.7827, gap 0.0000 (CERTIFIED)      at 0.7 s
Both find the SAME optimum -> CBC had the optimal design but could not PROVE it; HiGHS proves it in
<1 s. >170x faster to certification. On the real field this took the full K=1-10 sweep (both
objectives) from "minutes with open gaps" to "8 s, every row gap 0.000".

## Why HiGHS
- Much stronger branch-and-bound + presolve than CBC; parallelises the tree search.
- Open-source (MIT), free, no license -- keeps the pipeline fully reproducible with no proprietary
  dependency (unlike Gurobi/CPLEX, which are faster still but license-gated).
- python-mip does not support a HiGHS backend, so the model was ported to HiGHS's native Python API
  (highspy). CBC is retained in solver_benchmark.py as a cross-check: both engines returning the same
  optimum independently validates the model is correct, not just fast.

## Decision
HiGHS (highspy) is the canonical Tier-1 solver. Reported certificates are HiGHS gap 0.000. CBC kept
only for the correctness cross-check. If ever a problem exceeds HiGHS, academic Gurobi is the fallback.
