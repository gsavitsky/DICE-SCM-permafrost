"""
Coupled forward model: DICE 2023 economy, two-box climate and welfare with a choice of carbon cycle and permafrost.

  cfg = dict(carbon='dfair' | 'scm',          # DICE 2023 carbon cycle, or the process-based SCM
             perma=None | dict(m=1.0, abrupt='smooth' | 'none'),
             a2=dice2023.A2_STD, rartp=dice2023.RARTP_STD,   # damage coefficient, rate of time preference
             fbx={...}, KM=145.0,                # SCM feedback switches and fertilisation constant (optional)
             E_fn=f(year), F_nonco2=f(year) or float)   # optional prescribed CO2 emissions and non-CO2 forcing

Configurations used here: DICE 2023 (dfair), DICE 2023 + permafrost (dfair + perma), DICE-C (scm), DICE-CP
(scm + perma). Controls are the DICE 2023 optimal paths; RK4 with dt = 0.05 yr from 2020.
Permafrost carbon is added to the atmosphere (SCM) or to the emissions entering the DFAIR reservoirs.
"""
import math

import numpy as np

import dice2023 as D
import dfair
import permafrost as PF
import scm

NE = 6                      # economy/climate states: Tb1, Tb2, K, cumulative industrial CO2, U_acc, F_GHG


def init(cfg):
    eco = [0.1477, 1.099454, 295.0, 633.5, 0.0, 0.518]
    if cfg['carbon'] == 'dfair':
        carb = list(dfair.INIT_2020)
    else:
        L = scm.state_2020(cfg.get('fbx'), cfg.get('KM', scm.KM))
        carb = [L['Ca'], L['Cv'], L['Cs'], L['Dcum']] + [L[k] for k in scm.XKEYS]
    perm = (PF.state_2020(cfg['perma']) + [0.0]) if cfg.get('perma') else [0.0] * 5
    return np.array(eco + carb + perm, dtype=float)


def rhs(S, y, T200, cfg, pulse_gtc=0.0):
    Tb1, Tb2, K, Ctot, Uacc, FGHG = S[:NE]
    L, aL, sig, eland, Fmisc, FGHGexo, CO2EGHGabateB, pback, remiss, miu, s = D.exo(y)
    T = Tb1 + Tb2
    Yg = aL * (L / 1000.0) ** (1 - D.gamma_) * K ** D.gamma_
    Dfrac = cfg['a2'] * T ** D.a3_
    Cab = Yg * sig * remiss * pback / D.theta2_ / 1000.0 * miu ** D.theta2_
    Y = Yg * (1 - Dfrac) - Cab
    I = s * Y; Cons = Y - I
    ECo2 = (sig * Yg + eland) * (1 - miu)
    Vc = min(D.sd_cons_ ** 2 * (y - 2020), D.sd_cons_ ** 2 * 235)
    STP = 1.0 / (1 + cfg['rartp']) ** (y - 2020) / (1 + -0.5 * Vc * D.eta_ ** 2) ** (y - 2020)
    c = Cons * 1000.0 / L
    U = (c ** (1 - D.eta_) - 1) / (1 - D.eta_) - 1
    pi = len(S) - 5
    if cfg['carbon'] == 'dfair':
        R = S[NE:NE + 4]; Cum = S[NE + 4]
        MAT = dfair.Meq + R.sum()
    else:
        MAT = S[NE]
    Efn = cfg.get('E_fn')
    if 'F_nonco2' in cfg:
        Fnon = cfg['F_nonco2'](y) if callable(cfg['F_nonco2']) else cfg['F_nonco2']
    else:
        Fnon = Fmisc + FGHGexo
    F = D.fco22x_ * math.log(MAT / D.Meq_) / math.log(2) + Fnon
    dTb1 = (D.Teq1_ * F - Tb1) / D.d1_; dTb2 = (D.Teq2_ * F - Tb2) / D.d2_
    if cfg.get('perma'):
        dperm, Fmelt = PF.flows(T, dTb1 + dTb2, S[pi], S[pi + 1], S[pi + 2], S[pi + 3], T200, cfg['perma'])
    else:
        dperm, Fmelt = np.zeros(5), 0.0
    if cfg['carbon'] == 'dfair':
        E = (sum(np.atleast_1d(Efn(y))) if Efn else ECo2 / 3.667) + Fmelt + pulse_gtc
        a = dfair.alpha_of(Cum - (MAT - dfair.Meq), T)
        dcarb = np.concatenate([np.array(dfair.ESH) * E - R / (np.array(dfair.TAU) * a), [E]])
    else:
        Ca, Cv, Cs, Dc = S[NE:NE + 4]; X = S[NE + 4:NE + 16]
        e_ff = sig * Yg * (1 - miu) * 0.2727; e_luc = eland * (1 - miu) * 0.2727
        if Efn:                                  # prescribed emissions: PgC/yr, or (fossil, land use) PgC/yr
            _e = np.atleast_1d(Efn(y)); e_ff, e_luc = float(_e[0]), (float(_e[1]) if len(_e) > 1 else 0.0)
        Ph, Rp, Tu, Rs, dX, asg = scm.fluxes(Ca, Cv, Cs, Dc, X, T, cfg.get('fbx'), cfg.get('KM', scm.KM))
        dcarb = np.concatenate([[-Ph + Rp + Rs - asg + e_ff + e_luc + Fmelt + pulse_gtc,
                                 Ph - Rp - Tu - e_luc, Tu - Rs, e_luc], dX])
    deco = np.array([dTb1, dTb2, I - D.dK_ * K, ECo2 / 3.666, U * L * STP * D.scale1_,
                     D.Fc1_ * CO2EGHGabateB * (1 - miu) + math.log(D.Fc2_) * FGHG])
    return np.concatenate([deco, dcarb, dperm]), (c, STP, Y, Yg, Cab, Dfrac, MAT)


def run(cfg, dt=0.05, pulse=None, nyears=180):
    """Annual output 2020 .. 2020 + nyears. pulse = (year, GtCO2) released evenly through that year."""
    S = init(cfg); n = int(round(1 / dt)); lag = int(round(200 / dt))
    buf = [PF.WF * PF.T_hist(1820 + j * dt) for j in range(lag)]; csum = float(sum(buf)); t2m = 0.0
    out, cons = [], []
    pi = len(S) - 5
    for kk in range(nyears * n + 1):
        y = 2020 + kk * dt
        if kk % n == 0:
            _, diag = rhs(S, y, 0.0, cfg)
            out.append(S.copy()); cons.append(diag)
        if kk == nyears * n:
            break
        T200 = max(t2m, csum / lag); t2m = T200
        pg = (pulse[1] / 3.667) if (pulse and pulse[0] <= y < pulse[0] + 1) else 0.0
        k1 = rhs(S, y, T200, cfg, pg)[0]; k2 = rhs(S + dt / 2 * k1, y + dt / 2, T200, cfg, pg)[0]
        k3 = rhs(S + dt / 2 * k2, y + dt / 2, T200, cfg, pg)[0]; k4 = rhs(S + dt * k3, y + dt, T200, cfg, pg)[0]
        S = S + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        new = PF.WF * (S[0] + S[1]); csum += new - buf[kk]; buf.append(new)
    A = np.array(out); C = np.array(cons)
    d = dict(TATM=A[:, 0] + A[:, 1], Uacc=A[:, 4], release=A[:, pi + 4], cons=C[:, :2], thawed=A[:, pi + 1],
             Cperm=A[:, pi], Ctot=A[:, 3], Y=C[:, 2], Ygross=C[:, 3], Dfrac=C[:, 5], MAT=C[:, 6])
    if cfg['carbon'] == 'dfair':
        d['CO2ppm'] = (dfair.Meq + A[:, NE:NE + 4].sum(axis=1)) * scm.PPM
    else:
        d['CO2ppm'] = A[:, NE] * scm.PPM
        d['land'] = A[:, NE + 1] + A[:, NE + 2]; d['ocean'] = A[:, NE + 4:NE + 8].sum(axis=1) * 12e-15
    return d


def scc(cfg, year=2020, nyears=180, base=None):
    """Social cost of carbon (2019$ per tCO2) of a 1 GtCO2 pulse in `year`, welfare to 2020 + nyears:
    -1000 * dW/dE / (dW/dC), with dW/dC = 1000 * scale1 * STP * c^-eta averaged over the pulse year."""
    b = base or run(cfg, nyears=nyears)
    p = run(cfg, pulse=(year, 1.0), nyears=nyears)
    i = year - 2020
    mu = [1000 * D.scale1_ * b['cons'][j][1] * b['cons'][j][0] ** (-D.eta_) for j in (i, i + 1)]
    return -1000 * (p['Uacc'][-1] - b['Uacc'][-1]) / (0.5 * (mu[0] + mu[1]))
