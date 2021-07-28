#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.04.1

bin/nextflow \
  run . \
  -main-script src/predict_modality/workflows/run_baseline_benchmark/main.nf \
  -entry run_task1_benchmark \
  --publishDir output/ \
  -resume \
  --censor_dataset__max_mod1_columns 1000 \
  --censor_dataset__max_mod2_columns 1000
