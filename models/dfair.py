"""
DFAIR, the DICE 2023 carbon cycle (FaIR-based, Millar et al. 2017): four reservoirs (geological, deep ocean,
biosphere, ocean mixed layer) with fixed emission shares and decay timescales, scaled by a state-dependent factor
alpha such that the 100-yr integrated impulse response iIRF100 = r0 + irC * (cumulative uptake) + irT * T.
"""
import math

Meq = 588.0                                      # pre-industrial atmospheric carbon (PgC)
ESH = [0.2173, 0.224, 0.2824, 0.2763]            # emission shares: geological, deep ocean, biosphere, mixed layer
TAU = [1e6, 394.4, 36.53, 4.304]                 # decay timescales (yr)
PARAMS = dict(r0=32.4, irC=0.019, irT=4.165)     # iIRF100 terms; irC / irT can be rescaled for experiments


def alpha_of(cacc, T):
    """Timescale scaling alpha that gives iIRF100 = r0 + irC * cacc + irT * T (bisection)."""
    target = PARAMS['r0'] + PARAMS['irC'] * cacc + PARAMS['irT'] * T
    f = lambda a: sum(a * e * t * (1 - math.exp(-100 / (a * t))) for e, t in zip(ESH, TAU)) - target
    lo, hi = 0.1, 100.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


INIT_2020 = [150.093, 102.698, 39.534, 6.1865, 633.5]   # reservoirs R0-R3 (PgC above Meq), cumulative emissions
