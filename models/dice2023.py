"""
DICE 2023 (Barrage & Nordhaus 2024): economy, two-box climate and welfare parameters, exogenous drivers and the
DICE 2023 cost-benefit optimal control paths (emissions-control rate MIU and savings rate S).

Exogenous drivers (population, TFP, carbon intensity, land-use CO2, non-CO2 forcing) follow the DICE 2023 GAMS
formulas on the 5-year period grid and are linearly interpolated to annual values.
"""
import json
import math
import os

import numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "inputs")

P = dict(
    # climate (two-box) and forcing
    Teq1=0.324, Teq2=0.44, d1=236.0, d2=4.07, fco22x=3.93, Meq=588.0,
    # economy and welfare
    dK=0.0819, gamma=0.3, a3=2.0, theta2=2.6, eta=0.95, pi_=0.05,
    remiss0=1.4, remiss80=1.21, pback2050=515.0, scale1=0.00891061, sd_cons=0.01,
    Fc1=0.00955, Fc2=0.861,
)
gamma_, a3_, theta2_ = P['gamma'], P['a3'], P['theta2']
fco22x_, Meq_ = P['fco22x'], P['Meq']
Teq1_, Teq2_, d1_, d2_ = P['Teq1'], P['Teq2'], P['d1'], P['d2']
Fc1_, Fc2_ = P['Fc1'], P['Fc2']
dK_, sd_cons_, eta_, scale1_ = P['dK'], P['sd_cons'], P['eta'], P['scale1']

RARTP_STD = math.exp(0.001 + 0.5 * P['pi_']) - 1          # DICE 2023 risk-adjusted rate of time preference
A2_STD = 0.003467                                           # DICE 2023 quadratic damage coefficient


def _exogenous():
    gama, pop1, popadj, popasym = 0.300, 7752.9, 0.145, 10825.0
    q1, AL1, gA1, delA = 135.7, 5.84, 0.066, 0.0015
    gsigma1, delgsig, asymgsig, e1, miu1 = -0.015, 0.96, -0.005, 37.56, 0.05
    eland0, deland = 5.9, 0.1
    F_Misc2020, F_Misc2100, ECO2eGHGB2020, ECO2eGHGB2100 = -0.054, 0.265, 9.96, 15.5
    n = 81
    pv = np.arange(1, n + 1)
    L = np.zeros(n); L[0] = pop1
    for i in range(n - 1):
        L[i + 1] = L[i] * (popasym / L[i]) ** popadj
    gA = gA1 * np.exp(-delA * 5 * (pv - 1))
    aL = np.zeros(n); aL[0] = AL1
    for i in range(n - 1):
        aL[i + 1] = aL[i] / (1 - gA[i])
    gsig = np.minimum(gsigma1 * delgsig ** (pv - 1), asymgsig)
    sigma = np.zeros(n); sigma[0] = e1 / (q1 * (1 - miu1))
    for i in range(n - 1):
        sigma[i + 1] = sigma[i] * np.exp(5 * gsig[i])
    eland = eland0 * (1 - deland) ** (pv - 1)
    pw = lambda a, b: np.where(pv <= 16, a + ((b - a) / 16) * (pv - 1), b)
    return dict(years5=2020 + 5 * (pv - 1), L=L, aL=aL, sigma=sigma, eland=eland,
                F_Misc=pw(F_Misc2020, F_Misc2100), CO2E_GHGabateB=pw(ECO2eGHGB2020, ECO2eGHGB2100))


EXO = _exogenous()
_c = json.load(open(os.path.join(DATA, "dice2023_optimal_controls.json")))
_miu_x = np.linspace(_c["miu"]["xmin"], _c["miu"]["xmax"], len(_c["miu"]["ypts"])); _miu_y = np.array(_c["miu"]["ypts"])
_s_x = np.linspace(_c["s"]["xmin"], _c["s"]["xmax"], len(_c["s"]["ypts"])); _s_y = np.array(_c["s"]["ypts"])
_f = json.load(open(os.path.join(DATA, "dice2023_nonco2_fghg_forcing.json")))
_fg_x = np.linspace(_f["xmin"], _f["xmax"], len(_f["ypts"])); _fg_y = np.array(_f["ypts"])


def _interp(name, y):
    return float(np.interp(y, EXO['years5'], EXO[name]))


def _exo_at(y):
    pback = (P['pback2050'] * math.exp(-0.01 * (y - 2050)) if y - 2020 <= 30
             else P['pback2050'] * math.exp(-0.001 * (y - 2050)))
    remiss = P['remiss0'] + ((P['remiss80'] - P['remiss0']) / 80) * (y - 2020) if y <= 2100 else P['remiss80']
    return (_interp('L', y), _interp('aL', y), _interp('sigma', y), _interp('eland', y), _interp('F_Misc', y),
            float(np.interp(y, _fg_x, _fg_y)), _interp('CO2E_GHGabateB', y), pback, remiss,
            float(np.interp(y, _miu_x, _miu_y)), float(np.interp(y, _s_x, _s_y)))


_EXO_CACHE = {}


def exo(y):
    """(L, aL, sigma, eland, F_misc, F_GHG_exo, CO2e_GHG_abate_base, pbackstop, emissions ratio, MIU, S) in year y."""
    k = round(y, 6)
    v = _EXO_CACHE.get(k)
    if v is None:
        v = _EXO_CACHE[k] = _exo_at(y)
    return v
