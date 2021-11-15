#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.04.1


# bin/nextflow \
#   run . \
#   -main-script src/joint_embedding/workflows/censor_datasets/main.nf \
#   --datasets 'output/datasets/common/**.h5ad' \
#   --publishDir output/datasets/joint_embedding/ \
#   -resume


bin/nextflow \
  run . \
  -main-script src/joint_embedding/workflows/censor_datasets/main.nf \
  --datasets 'output/datasets_2021-11-08/phase1v2/common/**.h5ad' \
  --publishDir output/datasets_2021-11-08/phase1v2/joint_embedding/ \
  -resume \
  -c src/common/workflows/resource_labels_vhighmem.config \
  --censor_dataset__seed $SEED_SECRET \
  --censor_dataset__train_only false
bin/nextflow \
  run . \
  -main-script src/joint_embedding/workflows/censor_datasets/main.nf \
  --datasets 'output/datasets_2021-11-08/phase2_private/common/**.h5ad' \
  --publishDir output/datasets_2021-11-08/phase2/joint_embedding/ \
  -resume \
  -c src/common/workflows/resource_labels_vhighmem.config \
  --censor_dataset__seed $SEED_SECRET \
  --censor_dataset__train_only true
bin/nextflow \
  run . \
  -main-script src/joint_embedding/workflows/censor_datasets/main.nf \
  --datasets 'output/datasets_2021-11-08/phase2_private/common/**.h5ad' \
  --publishDir output/datasets_2021-11-08/phase2_private/joint_embedding/ \
  -resume \
  -c src/common/workflows/resource_labels_vhighmem.config \
  --censor_dataset__seed $SEED_SECRET \
  --censor_dataset__train_only false


