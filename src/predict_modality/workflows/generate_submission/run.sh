#!/bin/bash

set -ex

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# local repo & local datasets
export NXF_VER=21.04.1 

nextflow run \
  . \
  -main-script src/predict_modality/workflows/generate_submission/main.nf \
  --datasets 'output/public_datasets/predict_modality/**.h5ad' \
  --publishDir $HOME/Downloads/starter_kits/starter_kit-predict_modality-r/output/predictions/predict_modality \
  -resume