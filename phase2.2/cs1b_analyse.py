#!/usr/bin/env python3
"""
cs1b_analyse.py - S1.5: compares the first-order sphere-substrate energy from ORCA QM/MM single points with the
potential the A matrix uses (orca_vpot on the stored in vacuo density), at a contact, a mid and a far site.
  contact C: a from E(+-0.5), E(+-1) (Stage 1 runs s1a, s1e and S1.5 runs s1h, s1i), exact through q^3:
             D(q) = E(q) - E(-q) = 2 a q + 2 c q^3  ->  a = (8 D(0.5) - D(1)) / 6
  mid M, far F: a = (E(+1) - E(-1)) / 2 (exact through q^2; the q^3 term is negligible at 9 and 13 A)
  vpot: a_vpot = V(site) [Eh/e] x 627.5094740631
Criterion S1.5 (STAGE1B_CRITERIA.txt): |a_QMMM - a_vpot| <= 0.05 kcal/mol at all three sites.
Informative: the rigid model's electrostatics (AM1-BCC point charges of the committed force field, K sum q/r) at the
same sites, and the force-field-charge potential's share of the true one.
USAGE  python3 cs1b_analyse.py [phase2.2 folder on hpc1]
"""
import math, os, re, sys
P2 = sys.argv[1] if len(sys.argv) > 1 else "/home/18660916/system_development/phase2.2"
A, B = os.path.join(P2, "cs_stage1"), os.path.join(P2, "cs_stage1b")
H2KCAL = 627.5094740631
K = 332.0637133
def energy(path):
    t = open(path, errors="replace").read()
    if "ORCA TERMINATED NORMALLY" not in t:
        return None
    m = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", t)
    return float(m[-1]) if m else None
E = {"C+1": energy(os.path.join(A, "s1a_q1_N/job.out")), "C+0.5": energy(os.path.join(A, "s1e_q05_N/job.out")),
     "C-0.5": energy(os.path.join(B, "s1h_C_qm05/job.out")), "C-1": energy(os.path.join(B, "s1i_C_qm1/job.out")),
     "M+1": energy(os.path.join(B, "s1j_M_qp1/job.out")), "M-1": energy(os.path.join(B, "s1k_M_qm1/job.out")),
     "F+1": energy(os.path.join(B, "s1l_F_qp1/job.out")), "F-1": energy(os.path.join(B, "s1m_F_qm1/job.out"))}
rep = []
say = rep.append
say("S1.5 report (cs1b_analyse.py) - criterion from STAGE1B_CRITERIA.txt")
for k, v in E.items():
    say("  E(QM/MM) %-6s %s" % (k, "%.12f" % v if v is not None else "MISSING or not normal"))
if any(v is None for v in E.values()):
    say("S1.5 -> FAIL (a run is missing)"); open(os.path.join(B, "STAGE1B_REPORT.txt"), "w").write("\n".join(rep) + "\n"); print("\n".join(rep)); sys.exit(0)
D1 = (E["C+1"] - E["C-1"]) * H2KCAL
D5 = (E["C+0.5"] - E["C-0.5"]) * H2KCAL
aC = (8 * D5 - D1) / 6
cC = (D1 - 2 * aC) / 2
aM = (E["M+1"] - E["M-1"]) / 2 * H2KCAL
aF = (E["F+1"] - E["F-1"]) / 2 * H2KCAL
vl = [l.split() for l in open(os.path.join(B, "vpot_C_M_F.out")).read().splitlines() if l.strip()]
vals = [float(f[-1]) for f in vl if not (len(f) == 1 and f[0].isdigit())]
if len(vals) != 3:
    say("S1.5 -> FAIL (orca_vpot returned %d values, expected 3)" % len(vals))
    open(os.path.join(B, "STAGE1B_REPORT.txt"), "w").write("\n".join(rep) + "\n"); print("\n".join(rep)); sys.exit(0)
av = dict(zip("CMF", [v * H2KCAL for v in vals]))
# rigid-model (force-field charge) potential at the same sites
X = [[float(v) for v in l.split()[1:4]] for l in open(os.path.join(A, "s1a_q1_N/combined.xyz")).read().splitlines()[2:26]]
L = open(os.path.join(A, "s1a_q1_N/site.ORCAFF.prms")).read().splitlines()
i = L.index("$atoms")
q = [float(l.split()[2]) for l in L[i + 2:i + 2 + 24]]
S = {}
for l in open(os.path.join(B, "sites.tsv")).read().splitlines()[1:]:
    f = l.split("\t"); S[f[0]] = (float(f[1]), float(f[2]), float(f[3]), float(f[4]), int(f[5]))
aff = {k: K * sum(qj / math.dist(X[j], S[k][:3]) for j, qj in enumerate(q)) for k in "CMF"}
aq = {"C": aC, "M": aM, "F": aF}
say("")
say("site  nearest  beyond12  a_QMMM (kcal/mol/e)  a_vpot      difference   force-field charges (rigid model)")
ok = True
for k in "CMF":
    d = aq[k] - av[k]
    ok = ok and abs(d) <= 0.05
    say("  %s   %6.3f A   %2d/24    %12.4f      %12.4f   %+9.4f    %12.4f (%.1f%% of a_vpot)" % (
        k, S[k][3], S[k][4], aq[k], av[k], d, aff[k], 100 * aff[k] / av[k]))
say("contact-site cubic coefficient c = %.4f kcal/mol (the q^3 term removed by the symmetric fit)" % cC)
say("S1.5 QM/MM first-order sphere-substrate energy = A-matrix potential (orca_vpot, in vacuo density), within "
    "0.05 kcal/mol at contact, 9 A and 13 A -> %s" % ("PASS" if ok else "FAIL"))
if abs(aF) < 1.0 and abs(av["F"]) > 5.0:
    say("     note: the far-site QM/MM energy is near zero while the potential is not - a QM-MM cutoff is acting")
open(os.path.join(B, "STAGE1B_REPORT.txt"), "w").write("\n".join(rep) + "\n")
print("\n".join(rep))
