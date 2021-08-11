#!/bin/bash

# run prior to running this script:
# bin/viash_build -q common

set -ex

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

bin/viash run src/common/create_starter_kits/create_starter_kit/config.vsh.yaml -- \
  --input_dir src/predict_modality/starter_kit/starter_kit_r/ \
  --task predict_modality \
  --task_name "Predict Modality" \
  --language r \
  --language_name R \
  --block_starter 'par <- list(' \
  --evalai_phase 2276

bin/viash run src/common/create_starter_kits/create_starter_kit/config.vsh.yaml -- \
  --input_dir src/predict_modality/starter_kit/starter_kit_python/ \
  --task predict_modality \
  --task_name "Predict Modality" \
  --language python \
  --language_name Python \
  --block_starter 'par = dict(' \
  --evalai_phase 2276
  