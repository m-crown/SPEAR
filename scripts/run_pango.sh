#!/usr/bin/env bash -eu

CONDA_PATH=$(conda list -n base | awk 'NR == 1' | grep -oP '/\S+[^:]')
source ${CONDA_PATH}/etc/profile.d/conda.sh
conda activate pangolin
pangolin $1 --outdir $2
conda deactivate