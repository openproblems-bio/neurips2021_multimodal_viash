#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

bin/viash run src/common/datasets/dyngen/config.vsh.yaml -- \
  --id dyngen_bifurcating \
  --output resources/test/dyngen_bifurcating_antibody/dataset.h5ad \
  --plot resources/test/dyngen_bifurcating_antibody/plot.pdf \
  --backbone bifurcating_cycle \
  --num_cells 200 \
  --num_genes 150 \
  --num_simulations 30 \
  --num_threads 10 \
  --store_antibody

bin/viash run src/predict_modality/datasets/censor_task1/config.vsh.yaml -- \
  --input resources/test/dyngen_bifurcating_antibody/dataset.h5ad \
  --output resources/test/dyngen_bifurcating_antibody/dataset_task1_censor.h5ad
  
bin/viash run src/predict_modality/methods/baseline_randomforest/config.vsh.yaml -- \
  --input resources/test/dyngen_bifurcating_antibody/dataset_task1_censor.h5ad \
  --output resources/test/dyngen_bifurcating_antibody/dataset_task1_baseline.h5ad
