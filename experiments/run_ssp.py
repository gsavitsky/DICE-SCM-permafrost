"""
Emissions-driven SSPs, 2020-2100: RCMIP v5.1.0 CO2 emissions (MAGICC fossil & industrial, AFOLU) and non-CO2
effective radiative forcing (total minus CO2) into each carbon cycle, with and without permafrost, DICE 2023 two-box
climate. Reference: RCMIP SSP CO2 concentrations (MAGICC). Output: data/outputs/ssp_emissions_driven.csv
"""
import os
from multiprocessing import Pool
import pandas as pd
import _paths

SSPS = ['ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']
CONFIGS = {'DFAIR': ('dfair', False), 'DFAIR + permafrost': ('dfair', True), 'SCM': ('scm', False), 'SCM + permafrost': ('scm', True)}


def series(fn, var, scen):
    d = pd.read_csv(os.path.join(_paths.INPUTS, fn))
    r = d[(d.Region == 'World') & (d.Variable == var) & (d.Scenario == scen)]
    yrcols = [c for c in d.columns if c.isdigit()]
    s = r[yrcols].iloc[0].astype(float); s.index = s.index.astype(int)
    return s.interpolate()


def get(scen):
    ff = series('rcmip-emissions.csv', 'Emissions|CO2|MAGICC Fossil and Industrial', scen) / 3664.0
    lu = series('rcmip-emissions.csv', 'Emissions|CO2|MAGICC AFOLU', scen) / 3664.0
    erf = series('rcmip-radiative-forcing.csv', 'Effective Radiative Forcing', scen)
    eco2 = series('rcmip-radiative-forcing.csv', 'Effective Radiative Forcing|Anthropogenic|CO2', scen)
    conc = series('rcmip-concentrations.csv', 'Atmospheric Concentrations|CO2', scen)
    return ff, lu, erf - eco2, conc


def job(a):
    import model as M, dice2023 as D
    scen, name = a
    carbon, pf = CONFIGS[name]
    ff, lu, non, _ = get(scen)
    cfg = dict(carbon=carbon, a2=D.A2_STD, rartp=D.RARTP_STD, E_fn=lambda y: (float(ff[int(y)]), float(lu[int(y)])),
               F_nonco2=lambda y: float(non[int(y)]))
    if pf: cfg['perma'] = dict(m=1.0, abrupt='smooth')
    d = M.run(cfg, nyears=80)
    c = d['CO2ppm']
    return dict(model=name, scenario=scen, ppm2020=c[0], d2050=c[30] - c[0], d2100=c[80] - c[0], T2100=d['TATM'][80],
                release2100=d['release'][80])


if __name__ == "__main__":
    with Pool(8) as p:
        rows = p.map(job, [(s, n) for s in SSPS for n in CONFIGS], chunksize=1)
    for s in SSPS:
        conc = get(s)[3]
        rows.append(dict(model='RCMIP / MAGICC', scenario=s, ppm2020=conc[2020], d2050=conc[2050] - conc[2020],
                         d2100=conc[2100] - conc[2020]))
    D = pd.DataFrame(rows).sort_values(['scenario', 'model'])
    D.to_csv(os.path.join(_paths.OUTPUTS, "ssp_emissions_driven.csv"), index=False)
    print(D.pivot(index='model', columns='scenario', values='d2100').round(0).to_string())
