#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

out_dir=output/common_datasets/pbmc_1k_protein_v3
mkdir -p $out_dir

target/docker/common_datasets/download_10x_dataset/download_10x_dataset \
  --id pbmc_1k_protein_v3 \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5 \
  --output_rna $out_dir/pbmc_1k_protein_v3.normalise.output_rna.h5ad \
  --output_mod2 $out_dir/pbmc_1k_protein_v3.normalise.output_mod2.h5ad