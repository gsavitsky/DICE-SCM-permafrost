"""
Contribution of each carbon-cycle feedback to the SCC (Shapley decomposition) in DICE-CP.

Seven processes, each on or off (off = the process keeps its pre-industrial response, including over the 1800-2020
spin-up): ocean buffer saturation, ocean warming, CO2 fertilisation, photosynthesis temperature response, plant
respiration temperature response, soil respiration temperature response, permafrost. All 2^7 = 128 combinations are
run under the DICE 2023 optimal controls; SCC of 1 GtCO2 pulses in 2020 and 2080, welfare to 2200.
Shapley value of process i: phi_i = sum over subsets S of the other processes of |S|!(n-|S|-1)!/n! [v(S+i) - v(S)];
the seven values add up exactly to v(all on) - v(all off).
Outputs: data/outputs/shapley_runs.csv (all 128 runs), data/outputs/shapley_values.csv
"""
import itertools, math, os
from multiprocessing import Pool
import numpy as np, pandas as pd
import _paths

PROC = ['ocean', 'oceanT', 'fert', 'photoT', 'vegT', 'soilT', 'perma']
PULSES = (2020, 2080)


def job(c):
    import model as M, dice2023 as D
    cfg = dict(carbon='scm', a2=D.A2_STD, rartp=D.RARTP_STD, fbx=dict(zip(PROC[:6], c[:6])))
    if c[6]: cfg['perma'] = dict(m=1.0, abrupt='smooth')
    b = M.run(cfg)
    out = dict(zip(PROC, c), peak_T=b['TATM'].max(), T2100=b['TATM'][80], T2200=b['TATM'][180])
    for y in PULSES:
        out[f'SCC{y}'] = M.scc(cfg, y, base=b)
    return out


def shapley(v, n):
    phi = np.zeros(n)
    for i in range(n):
        others = [j for j in range(n) if j != i]
        for k in range(n):
            w = math.factorial(k) * math.factorial(n - k - 1) / math.factorial(n)
            for S in itertools.combinations(others, k):
                on = set(S)
                phi[i] += w * (v[tuple(1 if (j in on or j == i) else 0 for j in range(n))] -
                               v[tuple(1 if j in on else 0 for j in range(n))])
    return phi


if __name__ == "__main__":
    with Pool(8) as p:
        R = pd.DataFrame(p.map(job, list(itertools.product((0, 1), repeat=7)), chunksize=1))
    R.to_csv(os.path.join(_paths.OUTPUTS, "shapley_runs.csv"), index=False)
    rows = []
    for mtr in [f'SCC{y}' for y in PULSES] + ['peak_T', 'T2100', 'T2200']:
        v = {tuple(int(x) for x in r[PROC]): r[mtr] for _, r in R.iterrows()}
        phi = shapley(v, len(PROC))
        for p_, ph in zip(PROC, phi):
            rows.append(dict(metric=mtr, process=p_, shapley=ph, all_off=v[(0,) * 7], all_on=v[(1,) * 7]))
    S = pd.DataFrame(rows); S.to_csv(os.path.join(_paths.OUTPUTS, "shapley_values.csv"), index=False)
    print(S.pivot(index='process', columns='metric', values='shapley').loc[PROC].round(3).to_string())
