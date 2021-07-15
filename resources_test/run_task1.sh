#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

target_dir=target/docker
in_file=resources_test/common/pbmc_1k_protein_v3
out_file=resources_test/task1/pbmc_1k_protein_v3
mkdir -p `dirname $out_file`

$target_dir/predict_modality_datasets/prepare_task1_dataset/prepare_task1_dataset \
  --input_mod1 ${in_file}.output_rna.h5ad \
  --input_mod2 ${in_file}.output_mod2.h5ad \
  --output_mod1 ${out_file}.mod1.h5ad \
  --output_mod2 ${out_file}.mod2.h5ad \
  --output_solution ${out_file}.solution.h5ad
  
$target_dir/predict_modality_methods/baseline_randomforest/baseline_randomforest \
  --input_mod1 ${out_file}.mod1.h5ad \
  --input_mod2 ${out_file}.mod2.h5ad \
  --output ${out_file}.prediction.h5ad
  
$target_dir/predict_modality_metrics/calculate_task1_metrics/calculate_task1_metrics \
  --input_prediction ${out_file}.prediction.h5ad \
  --input_solution ${out_file}.solution.h5ad \
  --output ${out_file}.scores.h5ad
