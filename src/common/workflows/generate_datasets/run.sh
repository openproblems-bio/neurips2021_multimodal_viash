#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.04.1

bin/nextflow \
  run . \
  -main-script src/common/workflows/generate_datasets/main.nf \
  -entry generate_real_datasets \
  --publishDir output/public_datasets/common/ \
  -resume

bin/nextflow \
  run . \
  -main-script src/common/workflows/generate_datasets/main.nf \
  -entry generate_dyngen_datasets \
  --publishDir output/public_datasets/common/ \
  -resume
