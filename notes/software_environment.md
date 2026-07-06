# Software environment and module provenance

Records the actual software used, versions, and where modules were substituted
or worked around because of cluster runtime issues. Feeds the "Computational
details" part of the Methods chapter. Not a workflow step - reproducibility
metadata.

Cluster: hpc1.hpc (Stellenbosch HPC1). Scheduler: PBS.
Probed: 2026-07-06.

---

## Reference-paper software vs what this project used

The reference MD/QM-MM preparation protocol used **AMBER18**. AMBER18 is **not
available** on this cluster. The available AMBER is **AMBER22**, provided as two
modules:

| Module | Purpose | Tools provided |
|---|---|---|
| `app/amber22/22` | serial AMBER22 | antechamber, parmchk2, tleap, sander, pmemd |
| `app/ambermpi/22mpi` | MPI AMBER22 | pmemd.MPI, sander.MPI, pmemd, sander, cpptraj |

**Substitution rationale (for Methods):** AMBER22 replaces the paper's AMBER18.
The force-field and water-model choices were kept identical to the paper
(ff14SB protein, GAFF ligand, TIP3P water), so the version change affects the
MD engine implementation, not the underlying physics/parameters. Results are
therefore expected to be comparable. Do NOT claim AMBER18 was used.

`app/ambermpi/22mpi` loads requirements: `openmpi/4.1.1`, `compilers/gcc-9.4.0`.

---

## Tool availability (probed on login node)

Under `app/amber22/22` (serial) - loads `compilers/gcc-9.4.0`:

| Tool | Status |
|---|---|
| antechamber 22.0 | WORKS |
| parmchk2 | WORKS |
| tleap | WORKS |
| sander | login-node lib error (libquick.so) |
| pmemd | login-node lib error (libopenblas.so.0) |
| pmemd.MPI | not in serial module |
| cpptraj | BROKEN (see below) |

All structure-preparation and topology-building tools (antechamber, parmchk2,
tleap) work correctly. These cover the entire prep/parameterisation/solvation
phase.

---

## MD engines: login-node failure is NOT a real blocker

On the **login node**, the MD engines fail to load shared libraries:
- `pmemd` / `pmemd.MPI` / `sander.MPI`: `libopenblas.so.0: cannot open shared object file`
- `sander`: `libquick.so: failed to map segment`

**However, MD runs correctly on COMPUTE nodes via PBS.** Evidence:
- The accepted production segment was produced by
  `prod_seg001_1ns_pmemd_mpi32.pbs`, whose ONLY environment setup is
  `module purge; module load app/ambermpi/22mpi`, then
  `mpirun -np 32 pmemd.MPI ...`. No LD_LIBRARY_PATH edits, no conda, no extra
  library module - and it ran successfully.
- A full core-scaling benchmark tree (pmemd.MPI at 4/8/16/32 cores) also ran.

Conclusion: the missing `libopenblas.so.0` is a login-node-only gap; compute
nodes (allocated by PBS) carry the runtime. MD is never run interactively on the
login node, so this does not affect the workflow. No workaround was applied
(and none should be - the only system libopenblas.so.0 is an old mismatched copy
buried in an unrelated app, /apps/vrand/..., which must NOT be linked against
MD engines as it could silently corrupt numerics).

### Canonical MD run recipe (compute node, via PBS)
```
module purge
module load app/ambermpi/22mpi        # -> openmpi/4.1.1, gcc-9.4.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mpirun -np 32 pmemd.MPI -O -i in.mdin -o out.out -p top.prmtop \
       -c in.rst7 -r out.rst7 -x out.nc -inf out.mdinfo
```
PBS resource line used for accepted production:
`#PBS -l select=1:ncpus=32:mpiprocs=32:mem=64gb`, walltime 18:00:00.

---

## cpptraj: broken, separate issue (HPC-support item)

`cpptraj` (both serial and MPI modules) fails at runtime:
`libreadline.so.6: cannot open shared object file` and also `libquick.so`.
The system provides `libreadline.so.7` only (no `.so.6`).

- A project-local `libreadline.so.6 -> libreadline.so.7` symlink was tested in
  the original work and did NOT fix it (also needs libquick.so); the symlink was
  removed. Do not reintroduce this workaround.
- Whether cpptraj runs on compute nodes has not yet been confirmed the way pmemd
  has (via a successful job). TODO: test cpptraj in a PBS job, or request a
  working cpptraj from HPC support.
- Interim trajectory analysis: use Python/MDAnalysis/ParmEd-style tools or other
  available utilities instead of cpptraj until confirmed working.

---

## Local numpy note (non-cluster)

Local machine Python analysis used numpy 1.26.4 (in hpc1_home download). Not part
of the cluster MD/QM-MM runtime; listed only for completeness of local scripts.
