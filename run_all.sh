#!/usr/bin/env bash
# Run the whole lipidomics pilot analysis from the folder that holds the Excel file.
# Usage:  bash ~/lipidomics/run_all.sh [lipidomics-data.xlsx]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
X="${1:-lipidomics-data.xlsx}"
[ -f "$X" ] || { echo "Cannot find $X in $(pwd). cd to your data folder first."; exit 1; }
git -C "$HERE" pull -q origin claude/happy-fermat-kxnbea || echo "(could not pull; using existing code)"
run() { echo; echo "=== $1 ==="; shift; python "$HERE/$@"; }
run "1/6 Differential abundance (Mock/IFNB/IFNG)"         lipid_diff_analysis.py "$X" --outdir lipid_results
run "2/6 Class-group bar plots"                           lipid_group_plots.py   "$X" --outdir lipid_results/group_plots
run "3/6 Volcano-style plots (all + endolysosome)"        lipid_volcano.py       "$X" --outdir lipid_results/volcano
python "$HERE/lipid_volcano.py" "$X" --outdir lipid_results/volcano --classes LBPA GM1 GM2 SPH --tag endolysosome
run "4/6 Functional / process / localization sets"        lipid_set_analysis.py  "$X" --outdir lipid_results/sets
run "5/6 Lipid set enrichment + LION inputs (3 contrasts)" lipid_enrichment.py   "$X" --outdir lipid_results/enrichment
run "6/6 Indices, chain landscapes, FA routing, reaction map (3 contrasts)" lipid_indices.py "$X" --outdir lipid_results/indices
echo; echo "Done. Results in $(pwd)/lipid_results/"
