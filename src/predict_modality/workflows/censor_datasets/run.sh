#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.04.1

bin/nextflow \
  run . \
  -main-script src/predict_modality/workflows/censor_datasets/main.nf \
  --datasets 'output/public_datasets/common/**.h5ad' \
  --publishDir output/public_datasets/predict_modality/ \
  -resume

bin/nextflow \
  run . \
  -main-script src/predict_modality/workflows/censor_datasets/main.nf \
  --datasets 'output/private_datasets/common/**.h5ad' \
  --publishDir output/private_datasets/predict_modality/ \
  -resume
