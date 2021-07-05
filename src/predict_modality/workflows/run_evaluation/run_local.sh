#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"


export NXF_VER=21.04.1

# entry=dyngen_datasets
# tsv=src/common/datasets/dyngen/input.tsv
entry=public_10x_datasets
tsv=src/common/datasets/download_10x_dataset/input.tsv

bin/nextflow \
  run . \
  -main-script src/predict_modality/workflows/run_evaluation/main.nf \
  -entry $entry \
  --tsv $tsv \
  --publishDir output/ \
  -resume
