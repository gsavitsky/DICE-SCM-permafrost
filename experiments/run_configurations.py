"""
The four configurations (plus gradual-thaw-only permafrost) under the DICE 2023 optimal controls, 2020-2200:
peak warming, warming in 2100/2150/2200, permafrost carbon released, permafrost sensitivity (PgC released per degC of
warming by 2100) and the 2020 SCC (welfare to 2200). Output: data/outputs/configurations.csv
"""
import os
from multiprocessing import Pool
import pandas as pd
import _paths

CASES = {
    'DICE 2023': dict(carbon='dfair'),
    'DICE 2023 + permafrost (gradual only)': dict(carbon='dfair', perma=dict(m=1.0, abrupt='none')),
    'DICE 2023 + permafrost': dict(carbon='dfair', perma=dict(m=1.0, abrupt='smooth')),
    'DICE-C (SCM)': dict(carbon='scm'),
    'DICE-CP (SCM + permafrost, gradual only)': dict(carbon='scm', perma=dict(m=1.0, abrupt='none')),
    'DICE-CP (SCM + permafrost)': dict(carbon='scm', perma=dict(m=1.0, abrupt='smooth')),
}


def job(name):
    import model as M, dice2023 as D
    cfg = dict(CASES[name], a2=D.A2_STD, rartp=D.RARTP_STD)
    b = M.run(cfg)
    T = b['TATM']
    return dict(configuration=name, peak_T=T.max(), peak_year=2020 + int(T.argmax()), T2100=T[80], T2150=T[130],
                T2200=T[180], release_2100=b['release'][80], release_2200=b['release'][180],
                PgC_per_degC_2100=(b['release'][80] / T[80]) if cfg.get('perma') else float('nan'),
                SCC2020=M.scc(cfg, 2020, base=b))


if __name__ == "__main__":
    with Pool(6) as p:
        D = pd.DataFrame(p.map(job, list(CASES)))
    D.to_csv(os.path.join(_paths.OUTPUTS, "configurations.csv"), index=False)
    print(D.round(3).to_string(index=False))
