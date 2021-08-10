#!/bin/bash

set -ex 

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

target_dir=target/docker
out_file=resources_test/common/test_resource
mkdir -p `dirname $out_file`

$target_dir/common_datasets/download_totalvi_spleen_lymph/download_totalvi_spleen_lymph \
  --id test_resource \
  --input "https://github.com/YosefLab/totalVI_reproducibility/raw/master/data/spleen_lymph_111.h5ad" \
  --output_rna "${out_file}.tmp.output_rna.h5ad" \
  --output_mod2 "${out_file}.tmp.output_mod2.h5ad"

# stringent filtering to reduce the file size of test data
$target_dir/common_process_dataset/quality_control/quality_control \
  --input_rna "${out_file}.tmp.output_rna.h5ad" \
  --input_mod2 "${out_file}.tmp.output_mod2.h5ad" \
  --output_rna "${out_file}.output_rna.h5ad" \
  --output_mod2 "${out_file}.output_mod2.h5ad" \
  --min_counts_per_gene 10000 \
  --min_counts_per_cell 15000

rm "${out_file}.tmp.output_rna.h5ad" "${out_file}.tmp.output_mod2.h5ad"
