"""
spa_package_stats.py - the statistics of SP Amendment 2 (results 05_qmmm/24_sp_barriers/SP_AMENDMENT2.txt, A2.3), as
functions, used by spa_package.py. Python 3.6, standard library only (no numpy or scipy), deterministic.

  describe(values)            n, mean, sd (n - 1), sem, the 95% confidence interval of the mean (t-distribution, n - 1
                              degrees of freedom), median, min, max
  exp_average(values, T)      exponential average, -RT ln(mean exp(-E/RT)) (Cooper & Kastner 2014, their eq. 4;
                              SP_PLAN.txt), kcal/mol, computed stably
  cumulant(values, T)         Ryde 2017 eq. 5, mean - sd^2/(2RT), with the sample sd (n - 1) for sigma
  shapiro_wilk(values)        W and p (Royston 1995, Applied Statistics algorithm AS R94; n >= 3)
  skew_kurt(values)           adjusted Fisher-Pearson skewness G1 and excess kurtosis G2
  welch(a, b)                 mean(a) - mean(b) with its 95% confidence interval (Welch-Satterthwaite)
  ryde_table2(sigma_kj)       Ryde 2017 Table 2 (95% confidence): energies needed for theta = 10 and 20 kJ/mol, for the
                              exponential average and the cumulant approximation, at the first tabulated sigma at or
                              above the given one (None where the table has no entry)
  t_quantile(p, df), norm_cdf(z), norm_ppf(p)

Run "python3 spa_package_stats.py selftest" to compare every function with reference values computed independently
(scipy 1.18.1 on 18 October 2026; values embedded below).
"""
import math, sys

R_KCAL = 1.987204259e-3          # kcal/mol/K (CODATA R / 4184)
T_REF = 298.15                   # K, SP_PLAN.txt and SP_AMENDMENT2.txt A2.3
KJ = 4.184

def norm_cdf(z): return 0.5 * math.erfc(-z / math.sqrt(2.0))
def norm_ppf(p):
    """inverse of the standard normal cdf: Acklam's rational approximation refined by two Newton steps"""
    if not 0.0 < p < 1.0: raise ValueError("p must lie in (0, 1)")
    a = (-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02,
         -3.066479806614716e+01, 2.506628277459239e+00)
    b = (-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01,
         -1.328068155288572e+01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00, 2.938163982698783e+00)
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00)
    pl = 0.02425
    if p < pl:
        q = math.sqrt(-2 * math.log(p))
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    elif p > 1 - pl:
        q = math.sqrt(-2 * math.log(1 - p))
        x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    else:
        q = p - 0.5; r = q * q
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)
    for _ in range(2):                                   # Newton on the exact cdf
        x -= (norm_cdf(x) - p) / (math.exp(-0.5 * x * x) / math.sqrt(2 * math.pi))
    return x

def _betacf(a, b, x):
    """continued fraction of the incomplete beta function (modified Lentz)"""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny); h = d
    for m in range(1, 1000):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d; d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c; c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d; d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c; c = c if abs(c) > tiny else tiny
        de = d * c; h *= de
        if abs(de - 1.0) < 1e-15: return h
    raise ArithmeticError("incomplete beta did not converge")
def betainc(a, b, x):
    """regularised incomplete beta I_x(a, b)"""
    if x <= 0.0: return 0.0
    if x >= 1.0: return 1.0
    lbt = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1.0 - x)
    if x < (a + 1.0) / (a + b + 2.0): return math.exp(lbt) * _betacf(a, b, x) / a
    return 1.0 - math.exp(lbt) * _betacf(b, a, 1.0 - x) / b
def t_cdf(t, df):
    x = df / (df + t * t)
    tail = 0.5 * betainc(0.5 * df, 0.5, x)
    return 1.0 - tail if t >= 0 else tail
def t_quantile(p, df):
    """the t-distribution's p-quantile for p > 0.5 (bisection on the cdf to 1e-12)"""
    if not 0.5 < p < 1.0: raise ValueError("p must lie in (0.5, 1)")
    lo, hi = 0.0, 1.0
    while t_cdf(hi, df) < p: hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if t_cdf(mid, df) < p: lo = mid
        else: hi = mid
        if hi - lo < 1e-12: break
    return 0.5 * (lo + hi)

def describe(v):
    v = list(v); n = len(v)
    if n == 0: return None
    s = sorted(v); mean = sum(v) / n
    med = s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])
    out = {"n": n, "mean": mean, "median": med, "min": s[0], "max": s[-1]}
    if n >= 2:
        sd = math.sqrt(sum((x - mean) ** 2 for x in v) / (n - 1)); sem = sd / math.sqrt(n)
        t = t_quantile(0.975, n - 1)
        out.update(sd=sd, sem=sem, t95=t, ci95=(mean - t * sem, mean + t * sem))
    return out
def exp_average(v, T=T_REF):
    RT = R_KCAL * T; v = list(v)
    m = min(v)
    return m - RT * math.log(sum(math.exp(-(x - m) / RT) for x in v) / len(v))
def cumulant(v, T=T_REF):
    d = describe(v)
    if d is None or "sd" not in d: return None
    return d["mean"] - d["sd"] ** 2 / (2 * R_KCAL * T)
def skew_kurt(v):
    v = list(v); n = len(v)
    if n < 4: return None, None
    m = sum(v) / n
    m2 = sum((x - m) ** 2 for x in v) / n; m3 = sum((x - m) ** 3 for x in v) / n; m4 = sum((x - m) ** 4 for x in v) / n
    if m2 == 0: return None, None
    g1 = m3 / m2 ** 1.5; g2 = m4 / m2 ** 2 - 3.0
    G1 = g1 * math.sqrt(n * (n - 1)) / (n - 2)
    G2 = ((n + 1) * g2 + 6) * (n - 1) / ((n - 2) * (n - 3))
    return G1, G2
def _poly(c, x): return sum(ci * x ** i for i, ci in enumerate(c))
def shapiro_wilk(v):
    """W and p of the Shapiro-Wilk test, Royston (1995) AS R94 approximation, 3 <= n <= 5000"""
    x = sorted(v); n = len(x)
    if n < 3: return None, None
    if x[-1] - x[0] < 1e-12: return None, None
    mean = sum(x) / n; ssq = sum((xi - mean) ** 2 for xi in x)
    if n == 3:
        a = [-math.sqrt(0.5), 0.0, math.sqrt(0.5)]
    else:
        m = [norm_ppf((i - 0.375) / (n + 0.25)) for i in range(1, n + 1)]
        mtm = sum(mi * mi for mi in m); u = 1.0 / math.sqrt(n)
        an = m[-1] / math.sqrt(mtm) + _poly((0.0, 0.221157, -0.147981, -2.071190, 4.434685, -2.706056), u)
        if n > 5:
            an1 = m[-2] / math.sqrt(mtm) + _poly((0.0, 0.042981, -0.293762, -1.752461, 5.682633, -3.582633), u)
            eps = (mtm - 2 * m[-1] ** 2 - 2 * m[-2] ** 2) / (1 - 2 * an ** 2 - 2 * an1 ** 2)
            a = [m[i] / math.sqrt(eps) for i in range(n)]
            a[0], a[1], a[-2], a[-1] = -an, -an1, an1, an
        else:
            eps = (mtm - 2 * m[-1] ** 2) / (1 - 2 * an ** 2)
            a = [m[i] / math.sqrt(eps) for i in range(n)]
            a[0], a[-1] = -an, an
    W = sum(ai * xi for ai, xi in zip(a, x)) ** 2 / ssq
    W = min(W, 1.0)
    if W >= 1.0: return W, 1.0                          # a perfectly normal-looking sample: p = 1 (log(1 - W) undefined)
    if n == 3:
        p = max(0.0, min(1.0, 6.0 / math.pi * (math.asin(math.sqrt(W)) - math.asin(math.sqrt(0.75)))))
        return W, p
    if n <= 11:
        gamma = -2.273 + 0.459 * n
        mu = _poly((0.5440, -0.39978, 0.025054, -0.0006714), n)
        sigma = math.exp(_poly((1.3822, -0.77857, 0.062767, -0.0020322), n))
        if gamma - math.log(1 - W) <= 0: return W, 0.0
        z = (-math.log(gamma - math.log(1 - W)) - mu) / sigma
    else:
        ln = math.log(n)
        mu = _poly((-1.5861, -0.31082, -0.083751, 0.0038915), ln)
        sigma = math.exp(_poly((-0.4803, -0.082676, 0.0030302), ln))
        z = (math.log(1 - W) - mu) / sigma
    return W, 1.0 - norm_cdf(z)
def welch(a, b):
    da, db = describe(a), describe(b)
    va, vb = da["sd"] ** 2 / da["n"], db["sd"] ** 2 / db["n"]
    diff = da["mean"] - db["mean"]; se = math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (da["n"] - 1) + vb ** 2 / (db["n"] - 1))
    t = t_quantile(0.975, df)
    return {"diff": diff, "se": se, "df": df, "ci95": (diff - t * se, diff + t * se)}

# Ryde 2017, Table 2 (Just Accepted manuscript, PDF page 15, read from the rendered page): energies needed for the
# correct exponential average within theta with 95% confidence, as a function of sigma (kJ/mol).
# columns: sigma, EA theta=4, EA 10, EA 20, cumulant theta=4, cumulant 10, cumulant 20 (None: no entry in the table)
RYDE_T2 = (
    (1, 1, 1, 1, 2, 2, 2), (2, 2, 1, 1, 2, 2, 2), (3, 3, 1, 1, 5, 2, 2), (4, 8, 1, 1, 10, 3, 2), (5, 20, 2, 1, 18, 4, 2),
    (6, 42, 4, 1, 35, 6, 3), (7, 126, 8, 2, 55, 11, 4), (8, 372, 15, 3, 101, 15, 5), (9, 1560, 38, 4, 147, 24, 7),
    (10, 4180, 92, 7, 228, 35, 10), (11, 19900, 251, 12, 318, 50, 13), (12, 89300, 800, 28, 448, 68, 18),
    (13, 511000, 2660, 67, 593, 102, 24), (14, 3960000, 11800, 178, 809, 127, 32),
    (15, 28400000, 57200, 535, 1040, 167, 42), (16, None, 317000, 1780, 1330, 214, 53),
    (17, None, 2190000, 7410, 1680, 266, 66), (18, None, None, 33800, 2290, 340, 83),
    (19, None, None, 182000, 2650, 407, 106), (20, None, None, 1250000, 3170, 522, 127),
    (25, None, None, None, 7420, 1320, 292), (30, None, None, None, 15760, 2590, 627),
    (35, None, None, None, 30100, 4580, 1120), (40, None, None, None, 51000, 7670, 2070))
def ryde_table2(sigma_kj):
    for row in RYDE_T2:
        if row[0] >= sigma_kj:
            return {"sigma_row": row[0], "EA10": row[2], "EA20": row[3], "CA10": row[5], "CA20": row[6]}
    return {"sigma_row": None, "EA10": None, "EA20": None, "CA10": None, "CA20": None}

# reference values from scipy 1.18.1 (scipy.stats.t.ppf, shapiro, skew(bias=False), kurtosis(bias=False), norm.ppf,
# scipy.special.betainc), computed in the sandbox on 18 October 2026
#   (scipy.stats.t.interval / t.ppf / shapiro / skew(bias=False) / kurtosis(bias=False) / norm.ppf / ttest_ind(equal_var=
#   False).confidence_interval; the exponential average with scipy.special.logsumexp), 298.15 K
SELFTEST_DATA = {
    "A": [14.83, 12.95, 14.77, 10.08, 11.49, 13.96, 15.06, 14.94, 15.5, 14.55, 10.53, 12.21, 14.98, 16.29, 15.31, 13.31,
          15.26, 9.92, 16.04, 15.71, 16.46, 12.8, 15.71, 10.29, 9.88, 13.96, 17.57, 16.47, 10.3, 14.35],
    "B": [8.2, 10.91, 13.43, 8.4, 15.94, 10.69, 16.2, 13.26, 11.0, 6.21, 13.35, 13.09],
    "C": [12.05, 11.03, 10.43, 17.28, 11.07, 11.2, 13.38],
    "D": [14.03, 12.96, 16.88]}
REF = {
    'A_mean': 13.849333333333332, 'A_sd': 2.2796444570534926, 'A_ci_lo': 12.99810010344297,
    'A_ci_hi': 14.700566563223694, 'A_median': 14.66, 'A_ea': 11.061941176925787, 'A_ca': 9.463754550154057,
    'A_sw_W': 0.9086502743658328, 'A_sw_p': 0.013760782162092183, 'A_skew': -0.5624461146752607,
    'A_kurt': -0.8868903607569236,
    'B_mean': 11.723333333333331, 'B_sd': 3.0706301791276815, 'B_ci_lo': 9.77234799651295,
    'B_ci_hi': 13.674318670153712, 'B_median': 12.045, 'B_ea': 7.647294530010232, 'B_ca': 3.7663635404122786,
    'B_sw_W': 0.9499431283161931, 'B_sw_p': 0.6361533697055699, 'B_skew': -0.26340244838577787,
    'B_kurt': -0.5872622685248099,
    'C_mean': 12.348571428571429, 'C_sd': 2.3750609014747996, 'C_ci_lo': 10.152006451846038,
    'C_ci_hi': 14.54513640529682, 'C_median': 11.2, 'C_ea': 11.158410667139641, 'C_ca': 7.588185289294866,
    'C_sw_W': 0.774132618171533, 'C_sw_p': 0.02251712778429119, 'C_skew': 1.8816085756411245,
    'C_kurt': 3.5980613574189917,
    'D_mean': 14.623333333333335, 'D_sd': 2.0262362481540324, 'D_ci_lo': 9.589883456378402,
    'D_ci_hi': 19.656783210288268, 'D_median': 14.03, 'D_ea': 13.52009262791747, 'D_ca': 11.158575844915484,
    'D_sw_W': 0.9356899869285292, 'D_sw_p': 0.5103233249178379,
    'AB_welch_lo': 0.05097966249310382, 'AB_welch_hi': 4.201020337506899,
    't975_df2': 4.302652729749462, 't975_df9': 2.262157162798205, 't975_df19': 2.0930240544083087,
    't975_df29': 2.045229642132703, 't975_df12.7': 2.165567731368725,
    'nppf_1e-09': -5.9978070150076865, 'nppf_0.01': -2.3263478740408408, 'nppf_0.3': -0.5244005127080409,
    'nppf_0.975': 1.959963984540054}
def selftest(tol=1e-9):
    got = {}
    for k, v in SELFTEST_DATA.items():
        d = describe(v); W, p = shapiro_wilk(v)
        got.update({k + "_mean": d["mean"], k + "_sd": d["sd"], k + "_ci_lo": d["ci95"][0], k + "_ci_hi": d["ci95"][1],
                    k + "_median": d["median"], k + "_ea": exp_average(v), k + "_ca": cumulant(v),
                    k + "_sw_W": W, k + "_sw_p": p})
        if len(v) >= 4:
            g1, g2 = skew_kurt(v); got.update({k + "_skew": g1, k + "_kurt": g2})
    w = welch(SELFTEST_DATA["A"], SELFTEST_DATA["B"]); got.update(AB_welch_lo=w["ci95"][0], AB_welch_hi=w["ci95"][1])
    for df in (2, 9, 19, 29, 12.7): got["t975_df%g" % df] = t_quantile(0.975, df)
    for p in (1e-9, 0.01, 0.3, 0.975): got["nppf_%g" % p] = norm_ppf(p)
    bad = 0
    for k in sorted(REF):
        ok = k in got and abs(got[k] - REF[k]) <= tol
        bad += not ok
        print("%-14s %22.15g %22.15g  %s" % (k, got.get(k, float("nan")), REF[k], "ok" if ok else "DIFFERS"))
    if set(got) != set(REF): bad += 1; print("the computed and reference keys differ")
    print("selftest: %d values, %s (tolerance %g)" % (len(REF), "all agree" if bad == 0 else "%d differ" % bad, tol))
    return bad

if __name__ == "__main__":
    if sys.argv[1:] == ["selftest"]: sys.exit(1 if selftest() else 0)
    sys.exit(__doc__)
