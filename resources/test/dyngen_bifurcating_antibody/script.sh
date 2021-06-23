#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

OUTPUT_DIR=resources/test/dyngen_bifurcating_antibody

if [ ! -f $OUTPUT_DIR/dataset.h5ad ]; then
  bin/viash run src/common/datasets/dyngen/config.vsh.yaml -- \
    --id dyngen_bifurcating \
    --output $OUTPUT_DIR/dataset.h5ad \
    --plot $OUTPUT_DIR/plot.pdf \
    --backbone bifurcating_cycle \
    --num_cells 200 \
    --num_genes 150 \
    --num_simulations 30 \
    --num_threads 10 \
    --store_protein
fi

bin/viash run src/predict_modality/datasets/prepare_task1_dataset/config.vsh.yaml -- \
  --input $OUTPUT_DIR/dataset.h5ad \
  --output_censored $OUTPUT_DIR/dataset_task1_censored.h5ad \
  --output_solution $OUTPUT_DIR/dataset_task1_solution.h5ad

bin/viash run src/predict_modality/methods/baseline_randomforest/config.vsh.yaml -- \
  --input $OUTPUT_DIR/dataset_task1_censored.h5ad \
  --output $OUTPUT_DIR/dataset_task1_prediction_randomforest.h5ad
  
bin/viash run src/predict_modality/metrics/calculate_task1_metrics/config.vsh.yaml -- \
  --input_solution $OUTPUT_DIR/dataset_task1_solution.h5ad \
  --input_prediction $OUTPUT_DIR/dataset_task1_prediction_randomforest.h5ad \
  --output $OUTPUT_DIR/dataset_task1_prediction_randomforest_scores.h5ad
