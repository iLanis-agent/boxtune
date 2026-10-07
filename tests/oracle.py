#!/usr/bin/env python3
"""Independent oracle for BoxTune. Re-derives every published relation from scratch:
- sealed: exact -3dB of a 2nd-order high-pass from its magnitude polynomial (quadratic in u^2)
- vented: Small's 4th-order coefficients from h/alpha/QL, magnitude via cmath (independent
  arithmetic path from the JS manual re/im split), cutoff by brentq-style bisection on a fine grid
- port: Helmholtz resonator f = c/2pi*sqrt(A/(V*Leff)) inverted for L
- Keele: published constants evaluated with Python floats
Writes tests/expected.json consumed by tests/run_tests.js."""
import json, math, cmath, os

C = 343.0

def sealed_f3_ratio(qtc):
    # |H|^2 = u^4/((1-u^2)^2 + u^2/Q^2) = 1/2 -> y=u^2 solves y^2 + (2 - 1/Q^2)y - 1 = 0
    q2 = 1.0 / qtc**2
    disc = (2 - q2)**2 + 4
    return math.sqrt(((q2 - 2) + math.sqrt(disc)) / 2)

def sealed_from_qtc(fs, qts, vas, qtc):
    alpha = (qtc / qts)**2 - 1
    vb = vas / alpha
    fc = fs * qtc / qts
    return dict(Vb=vb, fc=fc, f3=fc * sealed_f3_ratio(qtc), alpha=alpha, Qtc=qtc)

def sealed_from_vb(fs, qts, vas, vb):
    alpha = vas / vb
    qtc = qts * math.sqrt(alpha + 1)
    fc = fs * math.sqrt(alpha + 1)
    return dict(Vb=vb, fc=fc, f3=fc * sealed_f3_ratio(qtc), alpha=alpha, Qtc=qtc)

def keele_optimum(fs, qts, vas):
    return dict(Vb=15 * vas * qts**2.87, Fb=0.42 * fs * qts**-0.9, F3=0.26 * fs * qts**-1.4)

def keele_given_vb(fs, qts, vas, vb):
    return dict(F3=math.sqrt(vas / vb) * fs, Fb=(vas / vb)**0.32 * fs)

def vented_G(u, fs, qts, vas, vb, fb, ql):
    """Small 4th-order high-pass, evaluated with complex arithmetic (independent path)."""
    h = fb / fs
    alpha = vas / vb
    sh = math.sqrt(h)
    a1 = (ql + h * qts) / (sh * ql * qts)
    a2 = (h + (alpha + 1 + h * h) * ql * qts) / (h * ql * qts)
    a3 = (h * ql + qts) / (sh * ql * qts)
    s = 1j * u
    return s**4 / (s**4 + a1 * s**3 + a2 * s**2 + a3 * s + 1)

def vented_analyze(fs, qts, vas, vb, fb, ql=7.0):
    def mag(u):
        return abs(vented_G(u, fs, qts, vas, vb, fb, ql))
    curve = []
    for i in range(240):
        f = 12 * (400 / 12) ** (i / 239)
        curve.append(dict(f=f, dB=20 * math.log10(mag(f / fb))))
    peak = 0.0; peakU = 1.0
    for j in range(601):
        u = 0.4 * (3 / 0.4) ** (j / 600)
        m = mag(u)
        if m > peak: peak, peakU = m, u
    def cut(target):
        prev_u, prev_m = 2.5, mag(2.5)
        for k in range(1, 4001):
            u = 2.5 * (0.02 / 2.5) ** (k / 4000)
            m = mag(u)
            if m < target <= prev_m:
                lo, hi = u, prev_u
                for _ in range(60):
                    mid = (lo + hi) / 2
                    if mag(mid) >= target: hi = mid
                    else: lo = mid
                return (lo + hi) / 2 * fb
            prev_u, prev_m = u, m
        return None
    samples=[curve[i] for i in range(0,240,20)]  # 12 grid samples keep expected.json small
    return dict(f3=cut(1 / math.sqrt(2)), f6=cut(0.5), f10=cut(10**-0.5),
                peakDb=20 * math.log10(peak), peakF=peakU * fb, samples=samples)

def port_length(vb_liters, fb, dia_mm, n=1):
    d = dia_mm / 1000.0
    a = math.pi * d * d / 4
    v = vb_liters / 1000.0
    leff = C**2 * n * a / (4 * math.pi**2 * fb**2 * v)
    L = leff - 0.732 * d
    return dict(lengthMm=L * 1000, leffMm=leff * 1000, endCorrMm=0.732 * d * 1000,
                areaCm2=n * a * 1e4)

DRIVERS = {
    "budget10": dict(fs=28.9, qts=0.38, vas=57.5),
    "compact65": dict(fs=45.0, qts=0.32, vas=14.0),
    "pa12": dict(fs=40.0, qts=0.24, vas=95.0),
    "b4textbook": dict(fs=30.0, qts=0.38268, vas=42.47),
}

out = {"sealed_table": [], "sealed_qtc": {}, "sealed_vb": {}, "keele": {}, "keele_vb": {},
       "vented": {}, "ports": {}}
# published Qtc -> f3/fc table (usenclosure / Dickason): engine must land within rounding
for qtc, pub in [(0.5, 1.55), (0.6, 1.21), (0.7, 1.0), (0.8, 0.9), (0.9, 0.83),
                 (1.0, 0.79), (1.1, 0.76), (1.2, 0.74), (1.3, 0.72), (1.4, 0.71)]:
    out["sealed_table"].append(dict(qtc=qtc, ratio=sealed_f3_ratio(qtc), published=pub))
for name, d in DRIVERS.items():
    out["sealed_qtc"][name] = sealed_from_qtc(d["fs"], d["qts"], d["vas"], 0.7071)
    out["sealed_vb"][name] = sealed_from_vb(d["fs"], d["qts"], d["vas"], 20.0)
    out["keele"][name] = keele_optimum(d["fs"], d["qts"], d["vas"])
    out["keele_vb"][name] = keele_given_vb(d["fs"], d["qts"], d["vas"], 30.0)
# vented cases: B4 golden (lossless), B4 with QL=7, keele box for budget10, boomy small box
out["vented"]["b4_lossless"] = vented_analyze(30, 0.38268, 42.47, 42.47/1.41421356, 30.0, 1e9)
out["vented"]["b4_ql7"] = vented_analyze(30, 0.38268, 42.47, 42.47/1.41421356, 30.0, 7.0)
out["vented"]["budget10_keele"] = vented_analyze(28.9, 0.38, 57.5, 53.7, 29.0, 7.0)
out["vented"]["boomy_small"] = vented_analyze(28.9, 0.38, 57.5, 18.0, 40.0, 7.0)
for spec in [(53.7, 29.0, 100, 1), (30.0, 35.0, 75, 1), (60.0, 28.0, 100, 2), (20.0, 40.0, 50, 1)]:
    out["ports"]["%gL_%gHz_%gmm_x%d" % spec] = port_length(*spec)

path = os.path.join(os.path.dirname(__file__), "expected.json")
with open(path, "w") as f:
    json.dump(out, f)
print("wrote", path)
print("b4_lossless f3/fs:", out["vented"]["b4_lossless"]["f3"] / 30)
print("b4_ql7 f3:", out["vented"]["b4_ql7"]["f3"], "peakDb:", out["vented"]["b4_ql7"]["peakDb"])
print("boomy peakDb:", out["vented"]["boomy_small"]["peakDb"], "at", out["vented"]["boomy_small"]["peakF"])
