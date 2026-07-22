# Known issues & failure modes register (address before finalising Tier-2)

Consolidated from Phase 1 close-out, the phase2b Tier-1 rebuild, and the Tier-2 method review.
Two kinds of entry live here: issues actually FACED in prior runs (sections D/E/F and B3/B4/C3 —
real, experienced), and ANTICIPATED Tier-2 failure modes from the method review that have not yet
bitten because Tier-2 has not run (most of section A, plus B1/B5/C1). Each entry: how it bit us OR
how it could bite Tier-2 / why it matters / the guard in place. Status: [OPEN] must be resolved
before the parallel batch; [ARGUED] a physical/documentary argument that must still be empirically
confirmed (treat as OPEN for the checklist); [CLOSED] guard in place, verify still holds; [WATCH]
no action now but monitor per-design.

---

## A. External-charge / field issues (the heart of Tier-2)

### A1. Geometric collapse / "implosion" under relaxation  [OPEN — critical]
Distinct from spill-out. Spill-out is a *static, diffuse-basis* effect (density detaching to a far
charge). Implosion is *dynamic*: during opt/NEB the NUCLEI move, and the atom carrying developing
negative charge (ether O3) is pulled toward the +1 site. The >=2 A safety margin is asserted only
at the INPUT geometry; nothing guarantees atoms stay >=2 A during relaxation. Risk scales with
field strength -> worst for max-lowering K=3,4 (strongest, closest-packed fields), not the gentle
K=1 distributed.
GUARD (to add): dynamic min(atom_i, charge_j) sentinel over ALL atoms x ALL charges, after every
endpoint opt AND along every NEB image; flag if any pair < ~2.0 A. HOMO-sign catches spill-out;
only this min-distance monitor catches implosion. Run Stage 0 on K=2 AND K=4 max-lowering (worst
case) before trusting the batch.
FALLBACK if it fires: native ECP-only dummy centre (Marefat Khah & Hattig 2020) at the charge
site supplies Pauli repulsion inside a plain ORCA job — no ASH/QM-MM switch needed.

### A2. Electron spill-out onto the +1 charge  [ARGUED — verify in Stage 0]
Physical argument (sound but not yet confirmed on THIS system): a non-diffuse basis (def2-SVP;
r2SCAN-3c's def2-mTZVPP likewise) physically cannot describe density detaching >=2 A beyond the
vdW surface, and CPCM eps=4 screens the bare Coulomb attraction and keeps the -2 bound. So bare
%pointcharges should be safe at +/-1 e / >=2 A. VERIFY empirically in Stage 0: HOMO stays negative
(bound dianion) + no Loewdin charge leaking toward a .pc site. Do NOT "fix" spill-out by adding
diffuse functions (makes it worse under an attractive external charge). Until Stage 0 passes this
is an argument, not a closed guard.

### A3. In-basis over-polarisation  [OPEN — add to Stage 0]
Separate from A1 and A2: even with a non-diffuse basis and no detachment, density can pile
unphysically onto the SUBSTRATE atoms nearest the +1. HOMO sign will NOT catch this.
GUARD: in Stage 0 also inspect Loewdin/Mulliken atomic charges on atoms near the +1 vs the bare
substrate; flag anomalous accumulation.

### A4. Walls cannot be bolted onto bare .pc charges in plain ORCA  [ARGUED — confirm once]
%LJCoefficients is a QM/MM-path keyword only (tied to the MM atom list / force-field files);
orca_pc is electrostatic-only; DoEQ restores charge-charge Coulomb, not a wall. So LJ walls need
full ORCA QM/MM or ASH. This is a strong inference from documentation STRUCTURE, NOT a verbatim
statement — the source research flagged it as such. CONFIRM with a one-off test (does
%LJCoefficients error or silently no-op in a bare %pointcharges job?) before relying on
"native = no walls possible". The ECP-dummy route (A1 fallback) sidesteps this entirely.

### A5. CPCM cavity must exclude the external charges  [CLOSED]
Tier-1 Dv was computed with charges ABSENT (cavity around the 24 solute atoms only). Keep the
Tier-2 cavity solute-only (ORCA builds it around real atoms, not .pc sites) so the field is not
dielectric-screened and Tier-1/Tier-2 stay consistent. Do NOT enclose charges in the cavity.

### A6. No DoEQ / no charge attenuation  [CLOSED]
Omit DoEQ: the charge-charge self-energy is a per-design constant (fixed charges) that cancels in
TS-minus-reactant and matches Tier-1 (which had no charges). Charge attenuation (Kalka/Stare
2025) is demoted — it rescales magnitudes and would corrupt the certified MILP field.

---

## B. Geometry, frame & reference-state issues

### B1. Fixed-field frame drift (rigid-body)  [OPEN — load-bearing]
The differential field assumes .pc sites stay at the exact MILP Cartesian coords relative to BOTH
endpoints in every step. Two ORCA behaviours break this:
 (i) Geometry opt loses the input orientation in cartesian->internal coordinate conversion.
     GUARD: %geom coordsys cartesian end (preserves the lab frame).
 (ii) NEB re-aligns images each step by quaternion superposition (Quatern, default on) and may
      re-centre (Fix_center) -> slides the solute relative to fixed charges.
     GUARD: %neb Quatern no end, no COM reset, no reorientation.
VERIFY the exact 6.0.1 keyword spelling (Quatern no / Fix_center) against the manual or a
one-image test — if mis-spelled ORCA may silently ignore it and corrupt the field on every design.
Plus a per-step endpoint frame-RMSD-vs-input sentinel (rigid-body drift). NOTE B1 is rigid-body
drift; A1 is a single atom migrating inward — DIFFERENT sentinels, need BOTH.

### B2. Endpoints must be re-optimised under each field  [CLOSED]
Reusing bare-substrate endpoints injects interpolation artefacts (the field shifts the minima) —
the GOCAT v1 lesson. Each design re-opts reactant AND product under its own field before path
generation. Rank by dE! = E(relaxed TS under field) - E(field-relaxed reactant under field),
per-design reference (not a common one — a common reference mixes in reactant-binding terms).

### B3. Mixed active-region barrier (Phase 1)  [CLOSED — do not regress]
Phase 1 briefly mixed a ~1959-atom reactant/product region with a 102-atom TS region -> an
apples-to-oranges barrier. Fixed by re-optimising R and P in the SAME 102-atom reduced region on
a COMMON frozen bulk (18c ts_start splice). Barrier +16.00 (supersedes the mixed +15.08).
Lesson for Tier-2: R, TS, P of a design must all be computed under the SAME field and SAME level.

### B4. Stale geometry propagation  [CLOSED — do not regress]
Old phase2 numbers were computed on the NEB-TS geometry (2.173/2.547) not the validated OptTS
(2.111/2.526); every *_realgeom result was silently wrong until audited. phase2b was rebuilt on the
validated consistent-region geometry. Lesson: never let energetics ride on a superseded geometry;
re-extract from committed converged structures.

### B5. Frame consistency across grid/reactant/product  [OPEN — one-time check]
Grid frame = reactant frame = product frame must hold (the splice should guarantee it). VERIFY
once before the batch: centroid + per-atom RMSD alignment of reactant vs product vs the frame the
grid/charges were built in. A silent offset corrupts every design's field.

### B6. CP2K chorismate is a placeholder  [WATCH]
The CP2K chorismate geometry is pedagogical; must be replaced with a properly optimised substrate
before any Phase-2 ENERGETIC result is treated as final. Tier-2 barriers inherit this caveat until
resolved.

---

## C. Path / TS issues

### C1. Scan / constraint indices are illustrative  [OPEN — trivial but silent]
Doc blocks show Scan B 2 3; the real breaking bond is O3-C4 = atoms 7-8 (0-based), forming bond
C1-C6 = 0-12. ORCA counts from 0. Copying verbatim yields a meaningless scan. Set and verify all
indices against the actual 24-atom geometry before submission.

### C2. Restrained single-bond scan is a diagnostic, not a scorer  [CLOSED]
A scan along one bond hard-codes a concerted assumption and can give discontinuous/hysteretic
barriers if the mechanism changes off-coordinate. Use it only to seed NEB and flag collapse; SCORE
with NEB-CI (+NEB-TS), which samples the whole path and detects mechanism change.

### C3. Spurious low-frequency imaginary modes  [WATCH]
Phase 1 saw ghost imaginary modes (~-18 cm-1) near the true TS mode; dismissed as numerical, not
extra transition vectors. In Tier-2, require exactly ONE imaginary mode on each located TS (freq
job) and treat tiny (~<50 cm-1) extras as ghosts, not mechanism.

---

## D. Solver / optimisation issues (Tier-1, carried as lessons)

### D1. CBC stalls on the distributed tail  [CLOSED]
python-mip/CBC returns SUBOPTIMAL with unclosed gaps on the high-K distributed problems and emits
"Coin0505I Presolved problem not optimal" spam. Two records: the earlier solver benchmark saw
severe stalls (K=7 gap 0.95; once 0.80 vs the true 0.77); THIS session's re-run saw milder stalls
(gaps 0.038 / 0.139 / 0.061 at K=6/7/9) — same conclusion either way. Swapped to HiGHS (highspy):
gap 0.000 on all K, both objectives, <1 s/solve. Certified gap is the whole differentiator vs
GOCAT — never report a stalled/SUBOPTIMAL row as certified.

### D2. Certified = proxy-optimal, NOT catalysis-optimal  [CLOSED — framing]
See DESIGN_DECISIONS.md. gap 0.000 proves optimality of the frozen-Dv PROXY. Three gaps to
catalysis: frozen density, frozen path (TS can move), over-optimisation. Tier-2 decides catalysis.

### D3. Max-lowering is monotone by construction  [CLOSED — framing]
"More charges -> lower barrier" is an artefact of the objective, not chemistry; barrier goes
unphysical (negative) past ~K=3-4. Max-lowering column is a diagnostic (ceiling + validity limit),
never a "how many charges" answer. Distributed load-spreading + Tier-2 ranking answer that.

### D4. Distributed designs cluster at the proxy floor  [WATCH — affects correlation]
Distributed K=1-4 all hold total ~-5.3 (band floor) -> they cluster at one proxy value and add
little spread to the Spearman correlation. Discriminating power comes from max-lowering K=1->4
spread + controls. Report the distributed cluster as a separate "equal-proxy -> equal-catalysis?"
sub-test, and compute rho on the K<=3 valid-regime subset separately from the full set.

---

## E. ORCA output parsing issues

### E1. Multiple energy lines per cycle  [CLOSED — do not regress]
QM/MM ORCA output prints ~5 energy lines/cycle. MUST grep FINAL SINGLE POINT ENERGY (QM/MM)
(~-1076) not (MM) (~-239). For bare-substrate Phase-2 single points the scale is different again
(~-836, QM-only). Confusing the scales silently corrupts every barrier. Always confirm the
magnitude matches the expected state (QM/MM ~-1076 enzyme; QM-only ~-836 bare substrate).

### E2. HOMO not labelled "HOMO"  [CLOSED]
ORCA prints ORBITAL ENERGIES with an OCC column; HOMO = last OCC>1.5 row, LUMO = first OCC~0.
A naive grep for "HOMO" finds nothing. Parse by occupancy transition (validated awk in 02/03).

### E3. orca_vpot container-member calling convention  [CLOSED]
orca_vpot takes: arg1 = .gbw, arg2 = the .scfp MEMBER name (relative, e.g. sp_reactant.scfp),
arg3 = points (Bohr), arg4 = output, arg5 = the .densities CONTAINER base (full path). arg2 is a
member name, arg5 the container — non-obvious; getting it wrong yields an empty/rubbish potential.

---

## F. HPC / environment & scaffolding issues

### F1. Login-node numpy segfault  [CLOSED]
numpy segfaults on the login node without thread-limit env vars. All HPC-run scripts are
pure-Python OR set OPENBLAS_NUM_THREADS=1 / OMP_NUM_THREADS=1 / MKL_NUM_THREADS=1 explicitly.

### F2. Workstation vs HPC path confusion  [CLOSED — recurring]
Running an HPC block on the workstation (or vice versa) fails with "FAIL missing ~/system_dev..."
because the paths only exist on one machine. Confirm the prompt (HPC: [18660916@hpc1]) before
pasting a machine-specific block. Geometries/grids must be scp'd to the HPC before HPC jobs.

### F3. Pre-3.8 Python on HPC  [CLOSED]
No math.dist() in HPC-bound scripts (manual sqrt of sum of squares); set thread env vars
explicitly. (Workstation venv is 3.8+ so math.dist is fine THERE only — do not copy such a script
to the HPC unchanged.)

### F4. Full exact-Hessian OptTS is memory-prohibitive  [CLOSED — Phase 1]
Full active-region exact Hessian ~221 GB > 126 GB node limit. Reduced-region substrate-only
Hybrid Hessian (144 displacements) is the validated method. Relevant if any Tier-2 TS needs a
Hessian at scale — keep the QM region small.

### F5. .gitignore blocks .npz / binary artefacts  [CLOSED]
grid_final.npz needed force-add past .gitignore. Check git status after committing binary results.

### F6. set -e ordering bug in scaffold scripts  [CLOSED]
An early phase2b scaffold aborted under set -e due to statement ordering (a check evaluated before
its inputs existed). Order guards AFTER the files/vars they test; re-tested before commit.

### F7. cpptraj cluster library issues  [WATCH]
cpptraj has cluster-wide library problems on this HPC (noted in Phase 1). If any Tier-2 step needs
cpptraj, expect linkage issues and prefer alternatives.

---

## G. Statistics / ranking integrity

### G1. Retain failures as penalty ranks, not exclusions  [CLOSED — design]
SCF fail / collapse / spill-out flagged designs stay in the ranking (penalty at worst end;
barrierless collapse at most-reactive end; average-rank ties). Excluding them = survivorship bias,
inflates rho. A proxy-rated-great design that then collapses is a proxy FAILURE and must count.

### G2. Small-n Spearman  [CLOSED — design]
n~12: critical |rho|=0.591 (p=0.05), powered only for rho>~0.7. Use EXACT-permutation p (not the
t-approx, inaccurate for n<30). Push toward n=15-19 with extra controls / K=5,6 max-lowering if
compute allows. Report rho (primary) + Pearson r + sign-agreement + scatter.

### G3. Level consistency Tier-1 vs Tier-2  [CLOSED — design]
FINAL ranking energies at the Tier-1 level (B3LYP-D3BJ/def2-SVP/CPCM eps=4). Changing level
between proxy and validation confounds rank changes (path relaxation vs method). r2SCAN-3c only
as a cheap geometry driver, B3LYP-D3BJ single point on top, gated by a rho>=0.9 calibration on
3-4 designs. GFN2-xTB excluded from ranking. r2SCAN-3c SIE risk for zwitterionic TS -> the
calibration gate is mandatory, do not skip.

---

## H. Provenance & data-hygiene issues (Phase 1, carried as lessons)

### H1. Provenance discrepancy in audit summaries  [CLOSED — lesson]
step07b_audit_summary.tsv once disagreed with its committed inputs; resolved by RECONSTRUCTING the
generating script and re-deriving. Lesson: committed scripts and their outputs must agree exactly;
if they don't, rebuild the generator, don't patch the output. Use neutral commit language when
correcting provenance.

### H2. AM1-BCC vs MOE template charges differ  [WATCH]
AM1-BCC charges differ from MOE template charges by up to 0.21 e on catalytically relevant atoms.
The choice must be disclosed and justified in Methods; do not silently mix charge schemes.

---

## Pre-batch checklist (the OPEN / ARGUED items, in order)
1. [A1/A3/A2] Stage 0 on K=2 AND K=4 max-lowering: HOMO sign (spill-out A2) + atomic-charge
   (over-polarisation A3) + dynamic min atom-charge distance (implosion A1). Gate the batch on all
   three passing.
2. [B1] Verify exact NEB frame-lock keyword (Quatern no / Fix_center) in ORCA 6.0.1; add
   coordsys cartesian; wire the frame-RMSD sentinel.
3. [B5] One-time frame-consistency check: grid = reactant = product frame.
4. [C1] Set + verify all scan/constraint atom indices (O3-C4 = 7-8, C1-C6 = 0-12; 0-based).
5. [A4] One-off confirm %LJCoefficients no-ops in a bare job (or just commit to ECP-dummy fallback).
6. Only then build the PBS array (Stage 1 calibration -> Stage 2 full batch -> Stage 3 scoring).
