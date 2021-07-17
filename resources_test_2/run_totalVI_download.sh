#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

echo "$REPO_ROOT"

target_dir=output
out_file=resources_test_2/common/
mkdir -p `dirname $out_file`

echo "spleen_lymph_111"
$target_dir/download_totalvi_spleen_lymph \
  --id spleen_lymph_111 \
  --input "https://github.com/YosefLab/totalVI_reproducibility/raw/master/data/spleen_lymph_111.h5ad" \
  --output_rna ${out_file}spleen_lymph_111_output_rna.h5ad \
  --output_mod2 ${out_file}spleen_lymph_111_output_mod2.h5ad

echo "spleen_lymph_206"
$target_dir/download_totalvi_spleen_lymph \
  --id spleen_lymph_206 \
  --input "https://github.com/YosefLab/totalVI_reproducibility/raw/master/data/spleen_lymph_206.h5ad" \
  --output_rna ${out_file}spleen_lymph_206_output_rna.h5ad \
  --output_mod2 ${out_file}spleen_lymph_206_output_mod2.h5ad

echo "totalVI_10x_malt_10k"
$target_dir/download_totalvi_10x \
  --id totalVI_10x_malt_10k \
  --input "https://github.com/YosefLab/totalVI_reproducibility/raw/master/data/malt_10k_protein_v3.h5ad" \
  --output_rna ${out_file}totalVI_10x_malt_10k_output_rna.h5ad \
  --output_mod2 ${out_file}totalVI_10x_malt_10k_output_mod2.h5ad
  
echo "totalVI_10x_pbmc_10k"
$target_dir/download_totalvi_10x \
  --id totalVI_10x_pbmc_10k \
  --input "https://github.com/YosefLab/totalVI_reproducibility/raw/master/data/pbmc_10k_protein_v3.h5ad" \
  --output_rna ${out_file}totalVI_10x_pbmc_10k_output_rna.h5ad \
  --output_mod2 ${out_file}totalVI_10x_pbmc_10k_output_mod2.h5ad

echo "totalVI_10x_pbmc_5k"
$target_dir/download_totalvi_10x \
  --id totalVI_10x_pbmc_5k \
  --input "https://github.com/YosefLab/totalVI_reproducibility/raw/master/data/pbmc_5k_protein_v3.h5ad" \
  --output_rna ${out_file}totalVI_10x_pbmc_5k_output_rna.h5ad \
  --output_mod2 ${out_file}totalVI_10x_pbmc_5k_output_mod2.h5ad
