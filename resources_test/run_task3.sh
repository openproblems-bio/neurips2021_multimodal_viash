#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

target_dir=target/docker
in_file=resources_test/common/pbmc_1k_protein_v3
out_file=resources_test/task1/pbmc_1k_protein_v3
mkdir -p `dirname $out_file`

$target_dir/match_modality_datasets/process_task3_dataset/process_task3_dataset \
  --input_rna ${in_file}.normalize.output_rna.h5ad \
  --input_mod2 ${in_file}.normalize.output_mod2.h5ad \
  --output_mod1 ${out_file}.mod1.h5ad \
  --output_mod2 ${out_file}.mod2.h5ad \
  --output_solution ${out_file}.solution.h5ad
  
$target_dir/match_modality_methods/baseline_pca/baseline_pca \
  --input_mod1 ${out_file}.mod1.h5ad \
  --input_mod2 ${out_file}.mod2.h5ad \
  --output ${out_file}.prediction.h5ad

