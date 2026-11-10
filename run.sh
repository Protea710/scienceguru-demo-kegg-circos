#!/usr/bin/env bash
# Rebuild the figure from the two public sources and check it.
# Usage: bash run.sh
set -euo pipefail
cd "$(dirname "$0")"

python3 code/01_prepare_data.py
python3 code/02_plot_circos.py
python3 code/03_verify.py
