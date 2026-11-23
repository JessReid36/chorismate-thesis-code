"""cs6_d52_expected.py - what Claude expects of DIAGNOSTIC 52: written by cs6_d52_expect_pyscf.py on 23 November 2026,
before D52's pre-registration. COMPUTED IN A SECOND PROGRAM (PySCF 2.15.1, numpy 2.5.3, scipy 1.18.1); not ORCA, not a
measurement and not part of the record. One electron, charges of +0.25 e at positive sites of B_n11_q0.25.
LIMIT: minus the lowest level in a large basis at the charges, eV. LIMIT2: the same in a still larger basis.
HF: minus the lowest level in set X2 or X3 at all six ghost atoms. CAM: (minus the energy with CAM-B3LYP in those
functions on a fine grid with exact integrals, the electrons that grid counts). CAMCHECK: one case on a finer grid.
The positions were read from d50/CAMTZX_B_n11_q0.25_X2A_di_R, sha256 of job.inp eb544c513a6ac5aa95617219331e795746c188f6aac6d0167e14d3157dac9c83
and of pc.pc 384f7bf8147c9c810053b8191ad9e8c1f8f47460713ce3ec3c69dfaa508434ef."""
SITE = {1124: (44.359719, 18.449173, 38.635519), 1073: (44.4024, 19.918764, 41.602409), 1107: (47.226943, 19.66848, 39.707342), 1359: (38.109868, 26.056989, 36.686294), 1005: (44.945806, 28.680872, 39.334093), 1052: (45.968418, 22.367003, 43.115654)}
GHOST = {1124: (44.359985, 18.439184, 38.635906), 1073: (44.402767, 19.913157, 41.610681), 1107: (47.234399, 19.662553, 39.710388), 1359: (38.103441, 26.064545, 36.685027), 1005: (44.946702, 28.690682, 39.33237), 1052: (45.975362, 22.364183, 43.122274)}
LIMIT = {"s1": 0.850355, "pc": 2.096745, "pf": 1.571243, "t3": 3.370500, "f4": 4.274655, "x6": 5.026069}
LIMIT2 = {"s1": 0.850356, "pc": 2.096747, "pf": 1.571246, "t3": 3.370503, "f4": 4.274658, "x6": 5.026073}
HF = {
    ("X2", "s1"): 0.847773, ("X2", "pc"): 2.096447, ("X2", "pf"): 1.570831, ("X2", "t3"): 3.369519, ("X2", "f4"): 4.273542, ("X2", "x6"): 5.024608,
    ("X3", "s1"): 0.850328, ("X3", "pc"): 2.096464, ("X3", "pf"): 1.571178, ("X3", "t3"): 3.369531, ("X3", "f4"): 4.273553, ("X3", "x6"): 5.024625,
}
CAM = {
    ("X2", "s1"): (0.899381, 1.000001), ("X2", "pc"): (2.160548, 1.000000), ("X2", "pf"): (1.675500, 1.000000), ("X2", "t3"): (3.432328, 1.000000), ("X2", "f4"): (4.341089, 1.000000), ("X2", "x6"): (5.091809, 1.000000),
    ("X3", "s1"): (0.903487, 1.000001), ("X3", "pc"): (2.160580, 1.000000), ("X3", "pf"): (1.675730, 1.000000), ("X3", "t3"): (3.432334, 1.000000), ("X3", "f4"): (4.341095, 1.000000), ("X3", "x6"): (5.091821, 1.000000),
}
CAMCHECK = {("X3", "f4"): (4.341095, 1.000000)}
