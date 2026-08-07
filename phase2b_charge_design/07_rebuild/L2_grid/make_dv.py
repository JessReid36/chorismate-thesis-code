#!/usr/bin/env python3
"""Layer 2 (part 2) -- dV = V_TS - V_R at every accepted grid site (proven phase2 dvpot recipe)."""
import os, subprocess
import numpy as np

HERE   = os.path.dirname(os.path.abspath(__file__))
INP    = os.path.join(HERE, "..", "inputs")
ORCA   = "/home/apps2/ORCA/6.0.1"
ORCA_BIN  = ORCA + "/orca"
ORCA_VPOT = ORCA + "/orca_vpot"
SIMPLE = "! B3LYP D3BJ def2-SVP def2/J RIJCOSX CPCM KeepDens"
BLOCKS = "%maxcore 3000\n%pal nprocs 8 end\n%cpcm epsilon 4.0 end\n%scf MaxIter 300 end"
ANG2BOHR = 1.8897259886


def read_grid():
    rows = [l.split("\t") for l in open(os.path.join(HERE, "grid_final.tsv")).read().splitlines()[1:]]
    xyz = np.array([[float(r[1]), float(r[2]), float(r[3])] for r in rows])
    return rows, xyz


def read_sub(geomfile):
    return open(os.path.join(INP, geomfile)).read().splitlines()[2:26]


def write_scf(base, sub):
    with open(os.path.join(HERE, base + ".inp"), "w") as f:
        f.write(SIMPLE + "\n" + BLOCKS + "\n* xyz -2 1\n")
        for l in sub:
            f.write(l.rstrip() + "\n")
        f.write("*\n")


def write_points_bohr(xyz_ang):
    with open(os.path.join(HERE, "points.xyz"), "w") as f:
        f.write("%d\n" % len(xyz_ang))
        for r in xyz_ang * ANG2BOHR:
            f.write("%.10f %.10f %.10f\n" % (r[0], r[1], r[2]))


def run_scf(base):
    with open(os.path.join(HERE, base + ".out"), "w") as fo:
        subprocess.run([ORCA_BIN, base + ".inp"], cwd=HERE, stdout=fo, stderr=subprocess.STDOUT, check=False)


def run_vpot(base):
    out = base + ".vout"
    subprocess.run([ORCA_VPOT, base + ".gbw", base + ".scfp", "points.xyz", out, base], cwd=HERE, check=False)
    return out


def parse_vpot(path, npts):
    vals = []
    for ln in open(os.path.join(HERE, path)).read().splitlines():
        p = ln.split()
        if len(p) >= 4:
            try:
                vals.append(float(p[-1]))
            except ValueError:
                pass
    assert len(vals) == npts, "parsed %d potentials, expected %d (check %s)" % (len(vals), npts, path)
    return np.array(vals)


def main():
    rows, xyz = read_grid()
    n = len(xyz)
    write_points_bohr(xyz)

    print("SCF + vpot: reactant ...")
    write_scf("dv_react", read_sub("reactant.xyz")); run_scf("dv_react")
    Vr = parse_vpot(run_vpot("dv_react"), n)
    print("SCF + vpot: TS ...")
    write_scf("dv_ts", read_sub("ts.xyz")); run_scf("dv_ts")
    Vts = parse_vpot(run_vpot("dv_ts"), n)

    dV = Vts - Vr
    with open(os.path.join(HERE, "grid_final.tsv"), "w") as f:
        f.write("idx\tx\ty\tz\tmin_dist\tnearest_el\tV_R\tV_TS\tDv\n")
        for r, vr, vts, dv in zip(rows, Vr, Vts, dV):
            f.write("\t".join(r[:6]) + "\t%.8f\t%.8f\t%.8f\n" % (vr, vts, dv))

    subR = np.array([[float(v) for v in l.split()[1:4]] for l in read_sub("reactant.xyz")])
    j = np.argmin(np.linalg.norm(xyz - subR[7], axis=1))
    print("dV: min %.5f max %.5f | |dV|max %.5f (a.u.), %d sites" % (dV.min(), dV.max(), abs(dV).max(), n))
    print("sign check: nearest site to ether-O (%.2f A) Dv=%+.6f -> +1 q*Dv=%+.6f %s"
          % (np.linalg.norm(xyz[j] - subR[7]), dV[j], dV[j],
             "STABILISING (correct)" if dV[j] < 0 else "destabilising (CHECK)"))
    print("rewrote grid_final.tsv with V_R, V_TS, Dv")


if __name__ == "__main__":
    main()
