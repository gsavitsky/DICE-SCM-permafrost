"""
Pairwise Shapley interaction indices (Grabisch & Roubens 1999) from the 128 runs in shapley_runs.csv:
  I_ij = sum over subsets S of the other processes of |S|!(n-|S|-2)!/(n-1)! [v(S+i+j) - v(S+i) - v(S+j) + v(S)]
plus, for each pair, the 2x2 contrast in each of the 32 settings of the other five processes: equal-weight mean,
value in the full model (others all on), range, and whether its sign is the same in all 32 settings.
Output: data/outputs/shapley_interactions.csv
"""
import itertools, math, os
import numpy as np, pandas as pd
import _paths

PROC = ['ocean', 'oceanT', 'fert', 'photoT', 'vegT', 'soilT', 'perma']
n = len(PROC)
R = pd.read_csv(os.path.join(_paths.OUTPUTS, "shapley_runs.csv"))
key = lambda on: tuple(1 if k in on else 0 for k in range(n))
rows = []
for mtr in [c for c in R.columns if c.startswith('SCC') or c in ('peak_T', 'T2100', 'T2200')]:
    v = {tuple(int(x) for x in r[PROC]): r[mtr] for _, r in R.iterrows()}
    for i, j in itertools.combinations(range(n), 2):
        oth = [k for k in range(n) if k not in (i, j)]
        tot, ds = 0.0, []
        for s in range(len(oth) + 1):
            w = math.factorial(s) * math.factorial(n - s - 2) / math.factorial(n - 1)
            for S in itertools.combinations(oth, s):
                S = set(S)
                tot += w * (v[key(S | {i, j})] - v[key(S | {i})] - v[key(S | {j})] + v[key(S)])
        for bg in itertools.product((0, 1), repeat=len(oth)):
            S = {o for o, x in zip(oth, bg) if x}
            ds.append(v[key(S | {i, j})] - v[key(S | {i})] - v[key(S | {j})] + v[key(S)])
        ds = np.array(ds)
        rows.append(dict(metric=mtr, i=PROC[i], j=PROC[j], interaction=tot, equal_weight=ds.mean(), others_all_on=ds[-1],
                         min=ds.min(), max=ds.max(), sign_consistent=bool((ds > 0).all() or (ds < 0).all())))
D = pd.DataFrame(rows); D.to_csv(os.path.join(_paths.OUTPUTS, "shapley_interactions.csv"), index=False)
print(D[D.metric == 'SCC2020'].assign(a=lambda t: t.interaction.abs()).sort_values('a', ascending=False)
      .head(8)[['i', 'j', 'interaction', 'min', 'max', 'sign_consistent']].round(2).to_string(index=False))
