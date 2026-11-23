"""
Permafrost carbon module following Hector v2.3pf (Woodard et al. 2021, GMD 14, 4751), with abrupt thaw following
Armstrong McKay et al. (2022).

Pools: frozen permafrost carbon (Cp), thawed permafrost carbon (Th); bookkeeping: cumulative abrupt (A) and gradual
(G) thaw, cumulative release to the atmosphere (R).
  frozen fraction  f(T) = 1 - Phi((ln(wf T) - mu) / sigma), wf = 2 (high-latitude amplification)
  gradual thaw     -df/dT * dT/dt * Cp (refreeze when cooling, limited by the thawed pool, timescale 2 yr)
  abrupt thaw      relaxes towards 0.5 * w(T) * G over 200 yr, w(T) = (1 + tanh((T - 1.5) / 0.5)) / 2
  respiration      J = m * 0.02 * (1 - f_static) * Th * Q10^(T200 / 10), T200 = 200-yr running mean of wf T,
                   never decreasing (as in Hector)
m is a strength multiplier on the respiration rate (m = 1 is Hector's value).
"""
import math

import numpy as np

WF, MU, SIG = 2.0, 1.67, 0.986
RATE, Q10, F_STATIC = 0.02, 2.0, 0.74
REFREEZE_TAU = 2.0
AB_T0, AB_W, AB_TAU, AB_FRAC = 1.5, 0.5, 200.0, 0.5
C0_PREINDUSTRIAL = 865.0
_SQ2PI = math.sqrt(2 * math.pi)


def dfrozen_dT(T):
    if T <= 1e-9:
        return 0.0
    z = (math.log(WF * T) - MU) / SIG
    return -math.exp(-0.5 * z * z) / _SQ2PI / (SIG * T)


def smax0(x, eps=1e-6):
    return 0.5 * (x + math.sqrt(x * x + eps))


def w_onset(T):
    return 0.5 * (1 + math.tanh((T - AB_T0) / AB_W))


def flows(T, dT, Cp, Th, A, G, T200, pm):
    """Tendencies [dCp, dTh, dA, dG, dR] (PgC/yr) and the flux to the atmosphere J.
    pm = dict(m=strength multiplier, abrupt='smooth'|'none')."""
    raw = -dfrozen_dT(T) * dT * Cp
    lim = -Th / REFREEZE_TAU
    grad = 0.5 * (raw + lim + math.sqrt((raw - lim) ** 2 + 1e-9))
    ab = smax0(AB_FRAC * w_onset(T) * smax0(G) - A) / AB_TAU if pm.get('abrupt') == 'smooth' else 0.0
    tot = grad + ab
    J = (1 - F_STATIC) * Th * RATE * pm.get('m', 1.0) * Q10 ** (T200 / 10.0)
    return np.array([-tot, tot - J, ab, grad, J]), J


def T_hist(y, k=0.02):
    return 0.0 if y < 1850 else 1.25 * (math.exp(k * (y - 1850)) - 1) / (math.exp(170 * k) - 1)


_CACHE = {}


def state_2020(pm):
    """[Cp, Th, A, G] in 2020 from a 1800-2020 spin-up from the pre-industrial pool (865 PgC, nothing thawed), forced
    by the historical warming reconstruction. Pre-2020 permafrost emissions are part of the observed historical carbon
    budget and are not added to the atmosphere."""
    key = (round(pm.get('m', 1.0), 6), pm.get('abrupt'))
    if key in _CACHE:
        return _CACHE[key]
    kk, dt = 0.02, 0.05
    dTdt = lambda y: 0.0 if y < 1850 else 1.25 * kk * math.exp(kk * (y - 1850)) / (math.exp(170 * kk) - 1)
    yrs = np.arange(1600, 2021)
    Ta = np.array([T_hist(y) for y in yrs])
    t200 = {int(y): WF * Ta[(yrs >= y - 200) & (yrs < y)].mean() for y in range(1800, 2021)}
    S = np.array([C0_PREINDUSTRIAL, 0.0, 0.0, 0.0])
    f = lambda S, y, T200: flows(T_hist(y), dTdt(y), S[0], S[1], S[2], S[3], T200, pm)[0][:4]
    for k in range(int(220 / dt)):
        y = 1800 + k * dt; T200 = t200[int(y)]
        k1 = f(S, y, T200); k2 = f(S + dt / 2 * k1, y + dt / 2, T200)
        k3 = f(S + dt / 2 * k2, y + dt / 2, T200); k4 = f(S + dt * k3, y + dt, T200)
        S = S + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    _CACHE[key] = [float(v) for v in S]
    return _CACHE[key]
