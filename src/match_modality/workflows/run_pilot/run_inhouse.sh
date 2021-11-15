#!/bin/bash

# run prior to running this script:
# bin/viash_build -q 'common'
# bin/viash_build -q 'match_modality' --max_threads 4

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.04.1

bin/nextflow \
  run . \
  -main-script src/match_modality/workflows/run_pilot/main.nf \
  -entry pilot_wf \
  --publishDir output/pilot/match_modality/ \
  -resume \
  -c src/common/workflows/resource_labels.config \
  --datasets 'output/datasets/match_modality/**.h5ad'

bin/nextflow \
  run . \
  -main-script src/match_modality/workflows/run_pilot/main.nf \
  -entry pilot_wf \
  --publishDir output/pilot/match_modality/ \
  -resume \
  -c src/common/workflows/resource_labels_highmem.config \
  --datasets 'output/datasets/match_modality/**.h5ad'