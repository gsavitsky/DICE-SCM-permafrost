---
output:
  pdf_document: default
  html_document: default
---
# Carbon-cycle feedbacks and permafrost in DICE 2023

This project asks **which carbon-cycle processes the social cost of carbon (SCC) depends on, and how much permafrost
thaw adds**. DICE 2023 represents the carbon cycle with DFAIR, an aggregate impulse-response emulator. Here DICE 2023
is coupled either to DFAIR or to a process-based simple carbon-cycle model (SCM), each with or without a permafrost
module. The SCM's processes can be switched to their pre-industrial response one at a time, which lets the SCC be
decomposed by process.

## Model

All configurations share the DICE 2023 economy, two-box temperature model, damages and welfare (Barrage & Nordhaus
2024) and are driven by the DICE 2023 cost-benefit optimal controls (emissions-control rate and savings rate).
Simulations run 2020-2200 (RK4, dt = 0.05 yr); welfare and the SCC are evaluated to 2200.

| Configuration | Carbon cycle | Permafrost |
|---|---|---|
| DICE 2023 | DFAIR | none |
| DICE 2023 + permafrost | DFAIR | yes |
| DICE-C | SCM | none |
| DICE-CP | SCM | yes |

- **DFAIR** (`models/dfair.py`): four reservoirs (geological, deep ocean, biosphere, ocean mixed layer) with fixed
  emission shares and timescales, scaled by a state-dependent factor so that the 100-yr integrated impulse response
  rises with cumulative uptake and warming (Millar et al. 2017).
- **SCM** (`models/scm.py`): atmosphere; vegetation and soil carbon (photosynthesis with CO2 fertilisation and a
  temperature response, Arrhenius plant respiration, Lloyd-Taylor soil respiration, land-use loss); a four-box ocean
  (warm and cold surface, intermediate, deep) with dissolved inorganic carbon, alkalinity and nitrate, full carbonate
  chemistry, a nutrient-driven biological pump and ocean-surface warming. Structure and parameters follow Lenton
  (2000, Tellus B) with Knox & McElroy (1984) ocean biology; K_M = 145. The 2020 state comes from a run starting at
  the pre-industrial steady state in 1800 under historical emissions (Global Carbon Budget from 1990) and warming.
- **Permafrost** (`models/permafrost.py`): Hector v2.3pf (Woodard et al. 2021): frozen and thawed pools, a
  lognormal frozen fraction of high-latitude warming, gradual thaw and refreeze, respiration of thawed carbon with a
  Q10 response, plus abrupt thaw following Armstrong McKay et al. (2022) (onset around 1.5 degC, 200-yr timescale,
  up to 50% additional thaw). Pools in 2020 come from an 1800-2020 spin-up from 865 PgC. Permafrost CO2 enters the
  atmosphere (SCM) or the DFAIR reservoirs.

### SCM feedback switches (`fbx` in `models/model.py`)
`ocean` (buffer saturation), `oceanT` (ocean warming), `fert` (CO2 fertilisation), `photoT`, `vegT`, `soilT`
(temperature responses of photosynthesis, plant and soil respiration), plus permafrost on/off. "Off" keeps the
process at its pre-industrial response, over the spin-up as well.

## Results so far (2020-2200, DICE 2023 optimal controls)

| | SCC 2020 ($/tCO2) | Peak warming (degC) | Warming 2100 / 2200 |
|---|---|---|---|
| DICE 2023 | 45.0 | 2.63 | 2.59 / 1.81 |
| DICE 2023 + permafrost | 48.9 | 2.81 | 2.73 / 2.21 |
| DICE-C | 45.6 | 2.59 | 2.56 / 1.53 |
| DICE-CP | 49.9 | 2.77 | 2.71 / 1.92 |

- **Evaluation.** The SCM's C4MIP feedbacks (1%/yr to 4xCO2) are land beta 0.81, ocean beta 0.86 PgC/ppm, land gamma
  -33, ocean gamma -21 PgC/degC, within the range of the 11 CMIP6 models in Arora et al. (2020). Splitting DFAIR by
  reservoir gives ocean beta 1.14 and ocean gamma -52, outside that range (fig. 1). Under emissions-driven SSPs, both
  carbon cycles with permafrost track RCMIP/MAGICC CO2 (e.g. SSP2-4.5: SCM 184, DFAIR 193, MAGICC 188 ppm rise
  2020-2100) and fall within AR6 assessed warming ranges (figs. 2-4).
- **Permafrost** releases about 105 PgC by 2100 on the DICE path (38 PgC per degC). This sits within the AR6 range
  (18, 3.1-41 PgC/degC) but toward its top, consistent with Hector's own results (figs. 5-6). Its contribution to
  warming after emissions stop (+0.10 degC at 50 yr) is consistent with MacDougall (2021) once his historical
  thaw-lag term is included (fig. 7).
- **SCC by process** (Shapley values, 2020 emission; fig. 8): with no carbon-cycle feedbacks the SCC is $42.8. Ocean
  buffer saturation +13.8, ocean warming +2.2, plant respiration +6.2, soil respiration +7.8 and permafrost +3.3
  raise it; CO2 fertilisation -14.2 and the photosynthesis temperature response -12.0 lower it; net +$7.1 ($49.9).
  Large opposing contributions nearly cancel. The same pattern holds for an emission in 2080 (fig. 9).
- **Interactions** (fig. 10): only the ocean-buffer interactions keep the same sign across all settings of the other
  processes. A saturating ocean amplifies the effect of every land process on the SCC, sinks and sources alike.

## Folder contents

| Path | Contents |
|---|---|
| `models/` | `dice2023.py` (economy, climate, exogenous drivers, optimal controls), `dfair.py`, `scm.py`, `permafrost.py`, `model.py` (coupled forward model, `run`, `scc`) |
| `experiments/` | scripts that produce every file in `data/outputs/`; `run_all.sh` runs them all and draws the figures |
| `data/inputs/` | DICE 2023 optimal controls and non-CO2 forcing, SCM pre-industrial steady state, Global Carbon Budget historical emissions, RCMIP v5.1.0 SSP data (world totals of the variables used), CMIP6 C4MIP values from Arora et al. (2020) Table A1 |
| `data/outputs/` | model results used by the figures (CSV) |
| `figures/` | `make_figures.R` and the ten figures |

| Output | Script | Figures |
|---|---|---|
| `configurations.csv` | `run_configurations.py` | 5, 6 |
| `c4mip_feedbacks.csv` | `run_c4mip.py` | 1 |
| `ssp_emissions_driven.csv` | `run_ssp.py` | 2, 3, 4, 6 |
| `zero_emissions_commitment.csv` | `run_zec.py` | 7 |
| `shapley_runs.csv`, `shapley_values.csv` | `run_shapley.py` (128 runs) | 8, 9 |
| `shapley_interactions.csv` | `shapley_interactions.py` | 10 |

## Running

Python 3 with numpy, scipy and pandas; R with readr, dplyr, tidyr and ggplot2.

```bash
cd experiments
./run_all.sh            # about 10 minutes on 8 cores; regenerates data/outputs and figures
```

A single run:
```python
import sys; sys.path.insert(0, "models")
import model, dice2023
cfg = dict(carbon="scm", perma=dict(m=1.0, abrupt="smooth"), a2=dice2023.A2_STD, rartp=dice2023.RARTP_STD)
d = model.run(cfg)                 # annual output 2020-2200: TATM, CO2ppm, release, Uacc, ...
print(d["TATM"].max(), model.scc(cfg, 2020))     # 2.768 degC, 49.88 $/tCO2
```
Switch processes with e.g. `cfg["fbx"] = dict(fert=0)`; scale permafrost respiration with `perma=dict(m=0.5, ...)`.

## Comparison sources used in the figures

Arora et al. (2020) Biogeosciences 17, 4173 (Table A1); Canadell et al. (2021) AR6 WG1 Ch. 5, Box 5.1; IPCC AR6 WG1
SPM Table SPM.1; Schädel et al. (2026) Commun. Earth Environ. 7, 306; Schuur et al. (2022) Annu. Rev. Environ.
Resour. 47, 343 (Table 1); Schneider von Deimling et al. (2015) and SROCC as assessed in AR6 Box 5.1; Woodard et al.
(2021) GMD 14, 4751 (Table 3); MacDougall (2021) Biogeosciences 18, 4937 (Table 2); RCMIP v5.1.0 (Nicholls et al.
2020; Meinshausen et al. 2020).

## Limits to keep in mind

- One damage function and discounting convention (DICE 2023 standard); the analysis is a forward simulation under
  the DICE 2023 optimal controls, not a re-optimisation.
- The SCM has fixed ocean circulation and biology; changes in stratification, overturning and the biological pump
  under warming are not represented.
- After 2200 the DICE 2023 control path's sustained net carbon removal pushes CO2 below pre-industrial in the SCM,
  so results are reported to 2200.
- Shares in the SCC decomposition depend on the SCM's calibration and on the permafrost and abrupt-thaw settings;
  the sweeps below are meant to test this.

## Suggested next analyses

Each sweep reruns the configurations (and, where useful, the 128-run Shapley decomposition) across a parameter range
anchored to published constraints. Report the SCC, peak and 2100/2200 warming, and the Shapley contributions.

### 1. CO2 fertilisation strength
- **What to vary:** the SCM's Michaelis-Menten constant K_M (`cfg["KM"]`; 145 at present). Lower K_M gives a
  stronger response to CO2.
- **Range:** choose K_M values that put the SCM's land beta (`experiments/run_c4mip.py`) across the CMIP6 spread,
  0.97 ± 0.40 PgC/ppm (Arora et al. 2020). Map K_M to land beta first, then pick, say, five points spanning the
  spread.
- **Questions:**
  - How do the SCC and the net carbon-cycle contribution change across the range?
  - Does the ranking of processes in the Shapley decomposition stay the same?
  - How does the fertilisation × ocean-buffer interaction scale?
- **Note:** the 2020 state is re-spun for each K_M (this is automatic in `scm.state_2020`). Check that 2020 CO2 and
  the historical land sink stay plausible.

### 2. Permafrost parameters
- **What to vary:**
  - respiration multiplier `m` (`perma=dict(m=...)`);
  - Q10 (2 at present);
  - static fraction f_static (0.74);
  - initial pool (865 PgC);
  - frozen-fraction parameters MU and SIG.
- **Constraints:** AR6 18 (3.1–41) PgC/°C (Canadell et al. 2021, Box 5.1); Schädel et al. (2026) gradual 28
  (10–36) PgC/°C; release by 2100 (Schneider von Deimling et al. 2015; SROCC); the ZEC contribution (MacDougall
  2021). m ≈ 0.5 gives about 22 PgC/°C, near the AR6 central value.
- **Design:** a one-at-a-time sweep first, then a small Latin-hypercube ensemble. Keep members consistent with the
  constraints above and report the SCC distribution.
- **Implementation:** Q10, f_static, MU, SIG and the pool size are module constants in `models/permafrost.py`. Move
  them into the `perma` dict, as `m` already is, so they can vary per run. Include them in the `state_2020` cache
  key so each member gets its own spin-up.

### 3. Abrupt-thaw functional form
The current form relaxes cumulative abrupt thaw towards 0.5 · w(T) · G over 200 yr, where
w(T) = (1 + tanh((T − 1.5)/0.5))/2. Abrupt thaw never reverses.
- **Parameters:**
  - onset temperature (1.0–2.5 °C, after Armstrong McKay et al. 2022);
  - transition width (0.25–1 °C);
  - timescale (50–300 yr);
  - maximum fraction of gradual thaw (0.2–0.5; or a fixed PgC per °C as in Schädel et al. 2026).
- **Shapes:** tanh (current), hard threshold (step), linear ramp, and an amount proportional to the remaining frozen
  pool instead of to gradual thaw.
- **Reversibility:** irreversible (current) or relaxing back below the threshold. This matters for the
  overshoot-and-return path in the DICE 2023 controls.
- **Metrics:** abrupt-thaw carbon by 2100/2200, the extra SCC, and the gap between gradual-only and gradual + abrupt
  results.
- **Implementation:** AB_T0, AB_W, AB_TAU and AB_FRAC in `models/permafrost.py` should move into the `perma` dict.
  Alternative shapes can be added as new values of `perma["abrupt"]` (currently `'smooth'` or `'none'`).
