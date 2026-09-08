# Candidate-grid build package

Builds two candidate-position grids from the same substrate geometries and
compares them:

- **`grid_poisson`** — Poisson-disk (blue-noise) elimination, the updated
  version of the existing pipeline
- **`grid_cvt_global`** and **`grid_cvt_per_shell`** — restricted centroidal
  Voronoi tessellation with Lloyd relaxation, the method your colleague uses

Both read the same signed distance field, extract the same shells, and sample
the same dense point cloud with the same seed. Only the point-selection step
differs, and the CVT grids are matched to the Poisson grid shell by shell, so
the comparison isolates the algorithm.

Everything is deterministic. Re-running with the same seed reproduces the same
grids exactly.

---

## What changed from the previous grid

The old pipeline used one parameter, `r_min = 1.5 Å`, for two unrelated jobs:
how finely a charge could be positioned, and how close two placed charges could
sit. Those are separated here.

**Resolution** is a property of the grid (`--r-min`, default 1.0 Å). Finer is
better, limited only by problem size.

**Separation** is a property of the species being placed — roughly 3.5 Å
centre-to-centre for formate, 4.5 Å for guanidinium — and belongs in the
optimiser as a constraint. These scripts do not decide it. They report how many
mutually compatible sites exist at each separation so the constraint can be set
knowingly.

Shells moved from 2/3/4 Å to **3/4/5 Å**. The offsets are measured from the
substrate van der Waals surface to a *point*; a molecular surrogate carries its
own radius, so a guanidinium carbon at 2 Å standoff puts its hydrogens through
the substrate surface. Pass `--shells 2.0 3.0 4.0` to reproduce the old
geometry for point-charge-only work.

Each site now also carries its measured standoff and the outward shell normal,
written to `*_sites.tsv`. The normal is what makes oriented surrogates
tractable: aligning a group's axis to the local normal removes two of the three
rotational degrees of freedom.

---

## Requirements

```bash
pip install numpy scipy scikit-image trimesh
```

On hpc1, inside the `ash` environment, or any Python 3.9+ with those four
packages. No ORCA and no electronic structure — this is pure geometry, and runs
in about two minutes on a laptop.

---

## Running it

Put these four scripts, `run_all.sh`, and the three substrate geometries in one
directory:

```bash
mkdir -p ~/gridbuild && cd ~/gridbuild
cp ~/Desktop/chorismate_thesis_code/phase2b_charge_design/01_geometry/reactant.xyz .
cp ~/Desktop/chorismate_thesis_code/phase2b_charge_design/01_geometry/ts.xyz .
cp ~/Desktop/chorismate_thesis_code/phase2b_charge_design/01_geometry/product.xyz .

# optional but recommended: calibrates the Dv-gradient column of the comparison
cp ~/Desktop/chorismate-thesis-results/phase2b_charge_design/03_dvpot/dv_grid.tsv .

bash run_all.sh 2>&1 | tee grid_build.log
```

### On the cluster instead

```bash
qsub run_all.pbs
```

Same scripts; the PBS wrapper only sets the thread environment numpy needs and
records the log. There is no advantage to running it on hpc1 unless you want
the record in the same place as everything else.

---

## What you get

| file | contents |
|---|---|
| `sdf_grid.npz` | union signed distance field over R + TS + P (72 atoms) |
| `grid_poisson.{npz,xyz}` | Poisson-disk grid |
| `grid_poisson_sites.tsv` | `idx, x, y, z, shell, standoff, nx, ny, nz` |
| `grid_cvt_global.*` | CVT, seeds partition the pooled cloud |
| `grid_cvt_per_shell.*` | CVT, each shell relaxed independently |
| `grid_comparison.txt` | the head-to-head table |

The `.xyz` files open directly in VMD, PyMOL or ChimeraX as dummy atoms, so you
can look at the two grids on the substrate side by side.

---

## Reading the comparison

**resolution (`res.mean`, `res.p95`, `res.max`)** — how far an arbitrary
position on the shells sits from the nearest available site. This is what limits
how well a charge can be positioned. Measured against a dense, independently
seeded reference set.

**Δv penalty** — the resolution error multiplied by the measured Δv gradient,
i.e. placement error expressed in kcal/mol. For scale, the best single site on
the old grid was −5.27 kcal/mol.

**min NN, pairs<0.8 Å, pairs<1.0 Å** — Poisson-disk guarantees a floor by
construction; Lloyd relaxation places no lower bound on any individual pair,
because it moves seeds to centroids and a centroid has no notion of a minimum
distance.

**CV** — the coefficient of variation of the spacing. This is the metric CVT
exists to minimise, so it is the fair test of that method's own claim.

**cond** — condition number of the inter-site 1/r matrix, the quadratic term in
the charge-selection problem. It depends on the smallest eigenvalue and is
therefore sensitive to a handful of near-degenerate pairs, so quote it alongside
the pair counts rather than alone.

**capacity** — how many mutually compatible sites exist at a given species
separation: the number of surrogates of that footprint that could coexist.

---

## Expected result

From a test run at 1448 points (365 / 479 / 604 on the 3, 4 and 5 Å shells):

| grid | N | min NN | CV | cond | pairs<1.0 Å | res.mean |
|---|---|---|---|---|---|---|
| Poisson | 1448 | 1.000 | 0.054 | 81,832 | 0 | 0.529 |
| CVT global | 1448 | 0.552 | 0.112 | 2,230,234 | 299 | 0.527 |
| CVT per-shell | 1448 | 0.552 | 0.112 | 255,473 | 300 | 0.527 |

The two methods cover the surface equally well — resolution is essentially
identical, 0.529 against 0.527 Å. The difference is entirely in the close pairs:
CVT produces 299 site pairs closer than 1.0 Å where Poisson produces none, and
worst-case inter-site coupling is 81% higher.

CVT also loses on CV, the metric it exists to optimise. The reason is that CVT
is a **surface** method while this grid is three concentric surfaces about 1 Å
apart — for many sites the nearest neighbour is on an adjacent shell, which
per-shell relaxation cannot see. Pooling the cloud (global mode) improves the
conditioning but cannot impose a floor.

**Caveat worth passing on:** if your colleague's grid is a *single* surface,
CVT is likely the better choice and none of this transfers. The result here is
specific to multiple closely spaced concentric shells.

---

## Tuning

```bash
# finer or coarser resolution
python3 02_grid_poisson.py --r-min 0.8 --shells 3.0 4.0 5.0 --density 10.0

# reproduce the old shell geometry for point-charge work
python3 02_grid_poisson.py --r-min 1.0 --shells 2.0 3.0 4.0 --density 6.0

# add a shell further out
python3 02_grid_poisson.py --r-min 1.0 --shells 3.0 4.0 5.0 6.0 --density 6.0
```

If you add a shell beyond about 8 Å, raise the margin in `01_build_sdf.py`;
the script will tell you if a requested shell falls outside the field.

Roughly, `--density` should be at least `2 / r_min²` so the dense cloud is not
the limiting factor: at `--r-min 1.0` use 6, at 0.8 use 10, at 0.6 use 18.

---

## After this

Δv must be re-evaluated on whichever grid you adopt — `orca_vpot` against the
existing reactant and TS densities, no new SCF. Every statistic in the current
Section X.4 refers to the old 331-point grid and will need replacing.

---

## Validation tests

Two additional checks, both pure geometry and both run from the same directory.

```bash
python3 05_convergence_test.py 2>&1 | tee convergence.txt
python3 06_envelope_test.py    2>&1 | tee envelope.txt
```

### 05 — voxel convergence

The signed distance field lives on a Cartesian lattice, and the lattice spacing
is a purely numerical parameter with no physical meaning. This halves it
(0.30 → 0.20 Å) and measures three things: shell area and offset, error against
the *analytic* surface, and the outcome for the finished grid.

The second measure is the decisive one. Because the union-of-spheres SDF is
analytic, the true shell is exactly SDF = d, so discretisation error can be
measured against ground truth rather than against the other mesh. Comparing the
two meshes to each other would be weaker and easy to get wrong — a nearest-vertex
query, for instance, reports the vertex spacing as disagreement even for two
identical surfaces.

Takes a few minutes; the fine lattice holds about 3.3× the voxels.

### 06 — envelope justification

Builds a grid from each single geometry and from the union, at identical
settings, then asks of each single-geometry grid how many of its sites fall
below the intended standoff once the omitted geometries are accounted for. Also
reports which geometry sets the envelope at each site of the union grid — if one
dominated, the union would be unnecessary.
