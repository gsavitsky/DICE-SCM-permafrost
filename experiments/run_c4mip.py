"""
C4MIP carbon-cycle feedback parameters, 1%/yr CO2 increase to 4xCO2 (140 yr) from pre-industrial:
beta (PgC/ppm) from the biogeochemically coupled run (no warming), gamma (PgC/degC) from the difference with the
coupled run (warming T = TCR log2(CO2/280), TCR = 1.8 degC), read at 4xCO2, as in Arora et al. (2020).
SCM: land and ocean carbon changes. DFAIR: compatible emissions from the prescribed concentration path; uptake by
reservoir = cumulative share of emissions - carbon still held; land = biosphere, ocean = deep ocean + mixed layer.
Output: data/outputs/c4mip_feedbacks.csv
"""
import math, os
import numpy as np, pandas as pd
import _paths
import scm, dfair

YEARS, TCR = 140, 1.8


def scm_run(coupled, dt=0.1):
    Ph, _, _, _ = scm.land(scm.CA_PI, 1.0, 1.0, 0.0, 0.0)
    P = scm.P
    Cv = Ph / (P['kr'] * P['kA'] * math.exp(-P['Ea'] / (P['R'] * 288)) + P['kt'])
    Cs = P['kt'] * Cv / (P['ksr'] * P['kB'] * math.exp(math.log(2.7) * (-(308.56 / (288 - 227.13)))))
    X = scm._YSS[1:].copy(); land0, ocean0 = Cv + Cs, X[0:4].sum() * 12e-15

    def f(st, t):
        ppm = 280 * 1.01 ** t; Ca = ppm / scm.PPM
        T = TCR * math.log2(ppm / 280) if coupled else 0.0
        Ph, Rp, Tu, Rs = scm.land(Ca, st[0], st[1], T, 0.0)
        dX, _ = scm.ocean_rhs(Ca, st[2:], T)
        return np.concatenate([[Ph - Rp - Tu, Tu - Rs], dX])
    st = np.concatenate([[Cv, Cs], X])
    for k in range(int(YEARS / dt)):
        t = k * dt
        k1 = f(st, t); k2 = f(st + dt / 2 * k1, t + dt / 2); k3 = f(st + dt / 2 * k2, t + dt / 2); k4 = f(st + dt * k3, t + dt)
        st = st + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return st[0] + st[1] - land0, st[2:6].sum() * 12e-15 - ocean0


def dfair_run(coupled, dt=0.05):
    R = np.zeros(4); cumE = 0.0; alloc = np.zeros(4)
    for k in range(int(YEARS / dt)):
        t = (k + 1) * dt
        ppm = 280 * 1.01 ** t; MATnew = dfair.Meq * ppm / 280.0
        T = TCR * math.log2(ppm / 280) if coupled else 0.0
        a = dfair.alpha_of(cumE - R.sum(), T)
        decay = R / (np.array(dfair.TAU) * a)
        E = ((MATnew - dfair.Meq) - (R - dt * decay).sum()) / (dt * sum(dfair.ESH))
        R = R + dt * (np.array(dfair.ESH) * E - decay); cumE += E * dt; alloc += np.array(dfair.ESH) * E * dt
    up = alloc - R
    return up[2], up[1] + up[3]                      # land (biosphere), ocean (deep + mixed layer)


if __name__ == "__main__":
    dppm = 280 * 1.01 ** YEARS - 280; dT = TCR * math.log2(1.01 ** YEARS)
    rows = []
    for name, fn in (("SCM", scm_run), ("DFAIR", dfair_run)):
        lb, ob = fn(False); lc, oc = fn(True)
        rows.append(dict(model=name, beta_land=lb / dppm, beta_ocean=ob / dppm, gamma_land=(lc - lb) / dT,
                         gamma_ocean=(oc - ob) / dT, beta_total=(lb + ob) / dppm, gamma_total=((lc + oc) - (lb + ob)) / dT))
    D = pd.DataFrame(rows); D.to_csv(os.path.join(_paths.OUTPUTS, "c4mip_feedbacks.csv"), index=False)
    print(D.round(3).to_string(index=False))
