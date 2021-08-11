#!/bin/bash

# run prior to running this script:
# bin/viash_build -q 'common'
# bin/viash_build -q 'joint_embedding' --max_threads 1

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.04.1

bin/nextflow \
  run . \
  -main-script src/joint_embedding/workflows/run_pilot/main.nf \
  -entry pilot_wf \
  --publishDir output/pilot/joint_embedding/ \
  -resume
