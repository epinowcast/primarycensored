# Reference log PMFs and gradients for the C++ uniform primary functions.
#
# Writes tests/testthat/testdata/uniform-lpmf-reference.csv: the uniform
# primary censored log PMF of [d - 1, d) for d = 1, ..., max_delay + 1, with
# pwindow = 1 and truncation to [0, max_delay + 1), and its gradient with
# respect to the unconstrained parameters theta used by the probe models in
# tests/testthat/helper-stan-cpp.R (mu and log(sigma) for lognormal, log of
# both parameters otherwise). Computed at 50 digits with mpmath.
#
# Usage: python data-raw/uniform-lpmf-reference.py
import csv
import mpmath as mp

mp.mp.dps = 50


def cdf_and_partial_expectation(dist, t, p1, p2):
    if dist == 1:
        mu, s = p1, p2
        z = (mp.log(t) - mu) / s
        return mp.ncdf(z), mp.exp(mu + s**2 / 2) * mp.ncdf(z - s)
    if dist == 2:
        k, b = p1, p2
        return (mp.gammainc(k, 0, b * t, regularized=True),
                k / b * mp.gammainc(k + 1, 0, b * t, regularized=True))
    k, lam = p1, p2
    x = (t / lam) ** k
    a = 1 + 1 / k
    return (1 - mp.exp(-x),
            lam * mp.gamma(a) * mp.gammainc(a, 0, x, regularized=True))


def G(dist, t, p1, p2):
    if t <= 0:
        return mp.mpf(0)
    F, M = cdf_and_partial_expectation(dist, t, p1, p2)
    return t * F - M


def lpmf(dist, max_delay, th1, th2):
    p1 = th1 if dist == 1 else mp.exp(th1)
    p2 = mp.exp(th2)
    g = [G(dist, t, p1, p2) for t in range(max_delay + 2)]
    cdf = [mp.mpf(0)] + [g[d] - g[d - 1] for d in range(1, max_delay + 2)]
    total = cdf[max_delay + 1]
    return [mp.log((cdf[d] - cdf[d - 1]) / total)
            for d in range(1, max_delay + 2)]


def num(x):
    """Format for R's read.csv(), which reads -Inf and NA but not -inf.

    Differences of log PMFs that are -inf are not numbers and are written
    as NA.
    """
    if isinstance(x, mp.mpc) or mp.isnan(x):
        return "NA"
    if mp.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    return mp.nstr(x, 20)


cases = [
    (1, "1.6", "0.42"), (1, "1.5", "0.5"), (1, "3.5", "0.05"),
    (2, "1.4", "0.38"), (2, "20", "4"), (2, "5", "1.2"),
    (3, "1.5", "5"), (3, "0.8", "3"), (3, "5", "0.8"),
]
h = mp.mpf("1e-20")
with open("tests/testthat/testdata/uniform-lpmf-reference.csv", "w",
          newline="") as f:
    out = csv.writer(f)
    out.writerow(["dist_id", "p1", "p2", "max_delay", "d", "lpmf",
                  "grad1", "grad2"])
    for dist, a, b in cases:
        th1 = mp.mpf(a) if dist == 1 else mp.log(mp.mpf(a))
        th2 = mp.log(mp.mpf(b))
        for max_delay in [14, 30]:
            v = lpmf(dist, max_delay, th1, th2)
            up1 = lpmf(dist, max_delay, th1 + h, th2)
            lo1 = lpmf(dist, max_delay, th1 - h, th2)
            up2 = lpmf(dist, max_delay, th1, th2 + h)
            lo2 = lpmf(dist, max_delay, th1, th2 - h)
            for i in range(max_delay + 1):
                out.writerow([
                    dist, a, b, max_delay, i + 1, num(v[i]),
                    num((up1[i] - lo1[i]) / (2 * h)),
                    num((up2[i] - lo2[i]) / (2 * h)),
                ])
