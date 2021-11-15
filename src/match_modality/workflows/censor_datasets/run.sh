#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.04.1



# bin/nextflow \
#   run . \
#   -main-script src/match_modality/workflows/censor_datasets/main.nf \
#   --datasets 'output/datasets/common/**.h5ad' \
#   --publishDir output/datasets/match_modality/ \
#   -resume \
#   -c src/common/workflows/resource_labels_highmem.config

bin/nextflow \
  run . \
  -main-script src/match_modality/workflows/censor_datasets/main.nf \
  --datasets 'output/datasets_2021-11-08/phase1v2/common/**.h5ad' \
  --publishDir output/datasets_2021-11-08/phase1v2/match_modality/ \
  -resume \
  -c src/common/workflows/resource_labels_highmem.config \
  --censor_dataset__seed $SEED_SECRET
bin/nextflow \
  run . \
  -main-script src/match_modality/workflows/censor_datasets/main.nf \
  --datasets 'output/datasets_2021-11-08/phase2_private/common/**.h5ad' \
  --publishDir output/datasets_2021-11-08/phase2_private/match_modality/ \
  -resume \
  -c src/common/workflows/resource_labels_highmem.config \
  --censor_dataset__seed $SEED_SECRET
# copy train files to phase2
find output/datasets_2021-11-08/phase2_private/match_modality/ -name '*_train_*' | sed 's#\(.*phase2\)_private\(.*\)\(/[^/]*\)$#mkdir -p \1\2; cp \1_private\2\3 \1\2\3#' | sh