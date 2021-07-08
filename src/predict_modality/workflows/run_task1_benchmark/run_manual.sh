#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

in_dir=output/common_datasets/pbmc_1k_protein_v3
out_dir=output/task1/pbmc_1k_protein_v3
mkdir -p $out_dir

target/docker/predict_modality_datasets/prepare_task1_dataset/prepare_task1_dataset \
  --input_rna $in_dir/pbmc_1k_protein_v3.normalize.output_rna.h5ad \
  --input_mod2 $in_dir/pbmc_1k_protein_v3.normalize.output_mod2.h5ad \
  --output_mod1 $out_dir/pbmc_1k_protein_v3.mod1.h5ad \
  --output_mod2 $out_dir/pbmc_1k_protein_v3.mod2.h5ad \
  --output_solution $out_dir/pbmc_1k_protein_v3.solution.h5ad