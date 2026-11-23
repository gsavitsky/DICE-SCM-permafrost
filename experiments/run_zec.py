"""
Permafrost contribution to the zero-emissions commitment: CO2 emissions grow 1%/yr from 11 PgC/yr in 2020 until
cumulative emissions reach 1000 PgC (DICE accounting, 633.5 PgC in 2020), then stop; non-CO2 forcing held at its 2020
value. Warming 50 and 100 yr after emissions stop, with and without permafrost, on each carbon cycle.
Output: data/outputs/zero_emissions_commitment.csv
"""
import math, os
from multiprocessing import Pool
import pandas as pd
import _paths


def job(a):
    import model as M, dice2023 as D
    carbon, pf = a
    cfg = dict(carbon=carbon, a2=D.A2_STD, rartp=D.RARTP_STD)
    if pf: cfg['perma'] = dict(m=1.0, abrupt='smooth')
    exo = D.exo(2020.0); F20 = exo[4] + exo[5]
    E20, g, cum0, total = 11.0, 0.01, 633.5, 1000.0
    tz = 2020 + math.log(1 + g * (total - cum0) / E20) / math.log(1 + g)
    d = M.run(dict(cfg, E_fn=lambda y: E20 * (1 + g) ** (y - 2020) if y < tz else 0.0, F_nonco2=F20), nyears=180)
    i0 = int(round(tz - 2020)); T = d['TATM']
    return dict(carbon='DFAIR' if carbon == 'dfair' else 'SCM', permafrost=pf, year_emissions_stop=tz,
                ZEC50=T[i0 + 50] - T[i0], ZEC100=T[i0 + 100] - T[i0])


if __name__ == "__main__":
    with Pool(4) as p:
        R = pd.DataFrame(p.map(job, [(c, pf) for c in ('dfair', 'scm') for pf in (False, True)]))
    rows = []
    for c in ('DFAIR', 'SCM'):
        a = R[(R.carbon == c) & ~R.permafrost].iloc[0]; b = R[(R.carbon == c) & R.permafrost].iloc[0]
        rows.append(dict(carbon=c, ZEC50_no_permafrost=a.ZEC50, ZEC50_with_permafrost=b.ZEC50,
                         permafrost_ZEC50=b.ZEC50 - a.ZEC50, ZEC100_no_permafrost=a.ZEC100,
                         ZEC100_with_permafrost=b.ZEC100, permafrost_ZEC100=b.ZEC100 - a.ZEC100))
    D = pd.DataFrame(rows); D.to_csv(os.path.join(_paths.OUTPUTS, "zero_emissions_commitment.csv"), index=False)
    print(D.round(3).to_string(index=False))
