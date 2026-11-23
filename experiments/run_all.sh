#!/bin/bash
# Regenerate every output in data/outputs (about 10 minutes on 8 cores), then the figures.
set -e
cd "$(dirname "$0")"
python3 run_configurations.py
python3 run_c4mip.py
python3 run_ssp.py
python3 run_zec.py
python3 run_shapley.py
python3 shapley_interactions.py
cd ../figures && Rscript make_figures.R
