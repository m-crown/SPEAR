#!/usr/bin/env bash
set -euo pipefail

CONDA_PATH=$(conda list -n base | awk 'NR == 1' | grep -oP '/\S+[^:]')
source ${CONDA_PATH}/etc/profile.d/conda.sh
conda activate pangolin
echo "Checking for Pangolin updates"
pangolin --update
echo "Running Pangolin"
pangolin ${1} --analysis-mode ${2} --outdir ${3} --max-ambig ${4} -t ${5} > ${6}
echo "Pangolin arguments: --analysis-mode ${2} --outdir ${3} --max-ambig ${4} -t ${5}" > ${3}/intermediate_output/pangolin_command.txt
pangolin --all-versions >> ${3}/intermediate_output/pangolin_command.txt
conda deactivate