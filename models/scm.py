"""
Simple carbon-cycle model (SCM) following Lenton (2000, Tellus B 52, 1159-1188), with ocean biology and carbonate
chemistry following Knox & McElroy (1984), version with ocean-surface warming (Lenton's Eq. 16) and K_M = 145.

  atmosphere (Ca); land: vegetation (Cv), soil (Cs), cumulative land-use loss (Dcum);
  ocean: four boxes (warm surface l, cold surface h, intermediate i, deep d), each with dissolved inorganic carbon,
  alkalinity and nitrate (mol). Air-sea exchange from surface pCO2 (carbonate system: Weiss 1974 K0, Goyet &
  Poisson 1989 K1/K2, Dickson 1990 KB, S = 34.7); circulation from Lenton's Fig. 1.

Feedback processes can be switched to their pre-industrial response (fbx dict; all on by default):
  ocean  : buffer saturation (off = surface pCO2 linear in DIC and alkalinity about the pre-industrial state)
  oceanT : ocean-surface warming (solubility and chemistry)
  fert   : CO2 fertilisation (off = photosynthesis sees pre-industrial CO2)
  photoT, vegT, soilT : temperature response of photosynthesis, plant respiration, soil respiration (off = no warming)
A switched-off process keeps its pre-industrial response through the 1800-2020 spin-up as well.
"""
import math
import os

import numpy as np
import pandas as pd
from scipy.optimize import brentq

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "inputs")

P = dict(rho=1027.0, Vl=3.06e16, Vh=5.4e15, Vi=3.6e17, Vd=9.7e17, kT=6.5e7, kU=4.87e7, kwi=1.25e7, kid=2e8,
         aO=3.6e14, Al=0.85, Ah=0.15, ka=1.773e20, fgtm=8.3259e13, ks=1752.0, Tl=295.15, Th=275.15,
         kp=0.184, kc=29.0, kr=0.092, kA=8.7039e9, Ea=54830.0, R=8.314, kt=0.092, ksr=0.0337, kB=157.072,
         kD=0.27, Cv0=550.0)
KM = 145.0
SY = 3.15576e7
S_ = 34.7
RHO = P['rho']
PPM = 1 / P['ka'] * P['fgtm'] * 1e6                       # ppm per PgC
V = dict(l=P['Vl'], h=P['Vh'], i=P['Vi'], d=P['Vd'])
AREA = dict(l=P['aO'] * P['Al'], h=P['aO'] * P['Ah'])
Tl0, Th0 = P['Tl'], P['Th']
kT, kU, kwi, kid = P['kT'], P['kU'], P['kwi'], P['kid']
B = 'lhid'
IDX = {b: k for k, b in enumerate(B)}
MASS = {b: V[b] * RHO for b in B}
_MBOX = np.array([MASS[b] for b in B])
ALPHA, CN, RAIN, F_ORG_I = 0.01, 7.0, 0.25, 0.3
Q = {('i', 'l'): kT + kwi, ('l', 'i'): kwi, ('l', 'h'): kT, ('i', 'h'): kU,
     ('h', 'd'): kT + kU, ('d', 'i'): kT + kU + kid, ('i', 'd'): kid}
XKEYS = [f"{l}_{b}" for l in ('DIC', 'TA', 'NO3') for b in B]


# ---------------- carbonate chemistry ----------------
def K0(T):
    return math.exp(-60.2409 + 93.4517 * (100 / T) + 23.3585 * math.log(T / 100)
                    + S_ * (0.023517 - 0.023656 * (T / 100) + 0.0047036 * (T / 100) ** 2))


def _K1K2(T):
    pK1 = 812.27 / T + 3.356 - 0.00171 * S_ * math.log(T) + 0.000091 * S_ ** 2
    pK2 = 1450.87 / T + 4.604 - 0.00385 * S_ * math.log(T) + 0.000182 * S_ ** 2
    return 10 ** -pK1, 10 ** -pK2


def _KB(T):
    s = S_
    return math.exp((-8966.90 - 2890.53 * s ** .5 - 77.942 * s + 1.728 * s ** 1.5 - 0.0996 * s ** 2) / T
                    + 148.0248 + 137.1942 * s ** .5 + 1.62142 * s
                    - (24.4344 + 25.085 * s ** .5 + 0.2474 * s) * math.log(T) + 0.053105 * s ** .5 * T)


TB = 0.000416 * S_ / 35.0


def pco2_surface(DIC, ALK, T):
    """DIC, ALK in mol/kg -> surface pCO2 in atm."""
    k1, k2 = _K1K2(T); kb = _KB(T)
    f = lambda H: DIC * (k1 * H + 2 * k1 * k2) / (H * H + k1 * H + k1 * k2) + TB * kb / (kb + H) - ALK
    H = brentq(f, 1e-11, 1e-5, xtol=1e-16)
    return DIC * H * H / (H * H + k1 * H + k1 * k2) / K0(T)


# ---------------- ocean ----------------
def _advect(c):
    d = np.zeros(4)
    for (src, dst), q in Q.items():
        f = q * SY * RHO * c[IDX[src]]
        d[IDX[src]] -= f; d[IDX[dst]] += f
    return d


def _biology(N):
    Nin_l = Q[('i', 'l')] * SY * RHO * N[2]
    Nin_h = Q[('l', 'h')] * SY * RHO * N[0] + Q[('i', 'h')] * SY * RHO * N[2]
    dD, dA, dN = np.zeros(4), np.zeros(4), np.zeros(4)
    for b, n in (('l', Nin_l), ('h', ALPHA * Nin_h)):
        corg = CN * n; caco3 = RAIN * corg; k = IDX[b]
        dD[k] -= corg + caco3; dA[k] += -2 * caco3 + n; dN[k] -= n
        dD[2] += F_ORG_I * corg; dN[2] += F_ORG_I * n; dA[2] -= F_ORG_I * n
        dD[3] += (1 - F_ORG_I) * corg + caco3; dN[3] += (1 - F_ORG_I) * n; dA[3] += 2 * caco3 - (1 - F_ORG_I) * n
    return dD, dA, dN


def ocean_rhs(Ca, X, dT=0.0):
    """X = [DIC(4), TA(4), NO3(4)] in mol. Returns dX (mol/yr) and net air-to-sea flux (PgC/yr)."""
    cD, cA, cN = X[0:4] / _MBOX, X[4:8] / _MBOX, X[8:12] / _MBOX
    pa = Ca / P['ka'] * P['fgtm']
    fl = {}
    for b, T in (('l', Tl0 + dT), ('h', Th0 + dT)):
        ps = pco2_surface(cD[IDX[b]], cA[IDX[b]], T)
        fl[b] = P['ks'] * AREA[b] * RHO * K0(T) * (pa - ps)
    dD, dA, dN = _biology(cN)
    dD = dD + _advect(cD); dA = dA + _advect(cA); dN = dN + _advect(cN)
    dD[0] += fl['l']; dD[1] += fl['h']
    return np.concatenate([dD, dA, dN]), (fl['l'] + fl['h']) * 12e-15


_YSS = np.load(os.path.join(DATA, "scm_preindustrial_steady_state.npy"))   # [Ca, DIC(4), TA(4), NO3(4)]
CA_PI = float(_YSS[0])
_C0 = {b: (_YSS[1 + IDX[b]] / _MBOX[IDX[b]], _YSS[5 + IDX[b]] / _MBOX[IDX[b]]) for b in 'lh'}
_LIN = {}
for _b, _T in (('l', Tl0), ('h', Th0)):
    _d, _a = _C0[_b]; _h = 1e-6
    _p0 = pco2_surface(_d, _a, _T)
    _LIN[_b] = (_p0, (pco2_surface(_d * (1 + _h), _a, _T) - pco2_surface(_d * (1 - _h), _a, _T)) / (2 * _d * _h),
                (pco2_surface(_d, _a * (1 + _h), _T) - pco2_surface(_d, _a * (1 - _h), _T)) / (2 * _a * _h))


def ocean_rhs_linear(Ca, X, dT=0.0):
    """ocean_rhs with surface pCO2 linear in DIC and alkalinity about the pre-industrial state (no buffer
    saturation); the temperature dependence of solubility and chemistry is kept."""
    cD, cA, cN = X[0:4] / _MBOX, X[4:8] / _MBOX, X[8:12] / _MBOX
    pa = Ca / P['ka'] * P['fgtm']
    fl = {}
    for b, T0 in (('l', Tl0), ('h', Th0)):
        k = IDX[b]; p0, gD, gA = _LIN[b]; T = T0 + dT
        ps = p0 + gD * (cD[k] - _C0[b][0]) + gA * (cA[k] - _C0[b][1])
        if dT:
            ps += pco2_surface(_C0[b][0], _C0[b][1], T) - p0
        fl[b] = P['ks'] * AREA[b] * RHO * K0(T) * (pa - ps)
    dD, dA, dN = _biology(cN)
    dD = dD + _advect(cD); dA = dA + _advect(cA); dN = dN + _advect(cN)
    dD[0] += fl['l']; dD[1] += fl['h']
    return np.concatenate([dD, dA, dN]), (fl['l'] + fl['h']) * 12e-15


# ---------------- land ----------------
def land(Ca, Cv, Cs, T, Dcum, km=KM):
    """Photosynthesis, plant respiration, litter turnover, soil respiration (PgC/yr); T = warming (degC)."""
    ppm = Ca * PPM
    Cvs = P['Cv0'] - P['kD'] * Dcum
    kMM = 1.478 * (251.0 / (120 + 251.0)) / (251.0 / (km + 251.0))      # pre-industrial photosynthesis fixed
    Photo = P['kp'] * Cvs * kMM * ((ppm - P['kc']) / (km + (ppm - P['kc']))) * ((15 + T) ** 2 * (25 - T) / 5625)
    PlantR = P['kr'] * Cv * P['kA'] * math.exp(-P['Ea'] / (P['R'] * (288 + T)))
    Turn = P['kt'] * Cv
    SoilR = P['ksr'] * Cs * P['kB'] * math.exp(math.log(2.7) * (-(308.56 / (288 + T - 227.13))))
    return Photo, PlantR, Turn, SoilR


def fluxes(Ca, Cv, Cs, Dc, X, T, fbx=None, km=KM):
    """Land fluxes and ocean tendency with feedback switches (see module docstring)."""
    fbx = fbx or {}
    on = lambda k: fbx.get(k, 1)
    Ph, _, _, _ = land(Ca if on('fert') else CA_PI, Cv, Cs, T if on('photoT') else 0.0, Dc, km)
    _, Rp, Tu, _ = land(Ca, Cv, Cs, T if on('vegT') else 0.0, Dc, km)
    _, _, _, Rs = land(Ca, Cv, Cs, T if on('soilT') else 0.0, Dc, km)
    dTo = T if on('oceanT') else 0.0
    dX, asg = ocean_rhs(Ca, X, dTo) if on('ocean') else ocean_rhs_linear(Ca, X, dTo)
    return Ph, Rp, Tu, Rs, dX, asg


# ---------------- 1800-2020 spin-up ----------------
def T_hist(y, k=0.02):
    """Historical warming reconstruction (degC above pre-industrial), 1.25 degC in 2020."""
    return 0.0 if y < 1850 else 1.25 * (math.exp(k * (y - 1850)) - 1) / (math.exp(170 * k) - 1)


# Lenton (2000) Table 3 historical emissions (total, land use; PgC/yr), used before 1990
_TOT = [(1800, .35), (1850, .4), (1947, 2.0), (1980, 6.7), (1990, 7.4)]
_LUC = [(1800, .35), (1850, .35), (1947, .8), (1980, 1.5), (1990, 1.4)]
_GCB = pd.read_csv(os.path.join(DATA, "historical_co2_emissions_gcb.csv"), index_col=0)   # MtCO2/yr


def _hist_emissions(y):
    """(fossil, land use) PgC/yr: Lenton Table 3 to 1990, Global Carbon Budget 1990-2019."""
    if y < 1990:
        xs, vs = zip(*_TOT); et = float(np.interp(y, xs, vs))
        xs, vs = zip(*_LUC); el = float(np.interp(y, xs, vs))
    else:
        el = _GCB.loc[y, 'land_use_change_co2'] / 3664.0; et = _GCB.loc[y, 'co2'] / 3664.0 + el
    return et - el, el


_STATE_CACHE = {}


def state_2020(fbx=None, km=KM):
    """SCM state in 2020: run from the pre-industrial steady state in 1800 under historical emissions and warming,
    with the same feedback switches and K_M as the run that uses it."""
    fbx = fbx or {}
    key = (tuple(sorted(fbx.items())), float(km))
    if key in _STATE_CACHE:
        return _STATE_CACHE[key]
    Ph, _, _, _ = land(_YSS[0], 1.0, 1.0, 0.0, 0.0, km)
    rP = P['kr'] * P['kA'] * math.exp(-P['Ea'] / (P['R'] * 288))
    Cv = Ph / (rP + P['kt'])
    Cs = P['kt'] * Cv / (P['ksr'] * P['kB'] * math.exp(math.log(2.7) * (-(308.56 / (288 - 227.13)))))
    S = np.concatenate([[_YSS[0], Cv, Cs, 0.0], _YSS[1:]]); dt = 0.1

    def f(S, T, eff, eluc):
        Ca, Cv, Cs, Dc = S[:4]
        Ph, Rp, Tu, Rs, dX, asg = fluxes(Ca, Cv, Cs, Dc, S[4:], T, fbx, km)
        return np.concatenate([[-Ph + Rp + Rs - asg + eff + eluc, Ph - Rp - Tu - eluc, Tu - Rs, eluc], dX])
    for y in range(1800, 2020):
        eff, el = _hist_emissions(y)
        for j in range(int(1 / dt)):
            a = (T_hist(y + j * dt), eff, el)
            k1 = f(S, *a); k2 = f(S + dt / 2 * k1, *a); k3 = f(S + dt / 2 * k2, *a); k4 = f(S + dt * k3, *a)
            S = S + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    st = dict(Ca=S[0], Cv=S[1], Cs=S[2], Dcum=S[3], **{k: S[4 + n] for n, k in enumerate(XKEYS)})
    _STATE_CACHE[key] = {k: float(v) for k, v in st.items()}
    return _STATE_CACHE[key]
