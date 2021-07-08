#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

target_dir=target/docker
out_file=resources_test/common/pbmc_1k_protein_v3
mkdir -p `dirname $out_file`

$target_dir/common_datasets/download_10x_dataset/download_10x_dataset \
  --id pbmc_1k_protein_v3 \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5 \
  --output_rna ${out_file}.output_rna.h5ad \
  --output_mod2 ${out_file}.output_mod2.h5ad

$target_dir/common_normalize/normalize/normalize \
  --input_rna ${out_file}.output_rna.h5ad \
  --input_mod2 ${out_file}.output_mod2.h5ad \
  --output_rna ${out_file}.normalize.output_rna.h5ad \
  --output_mod2 ${out_file}.normalize.output_mod2.h5ad \
  --min_counts_per_gene 1000 \
  --min_counts_per_cell 1000

rm ${out_file}.output_rna.h5ad ${out_file}.output_mod2.h5ad
