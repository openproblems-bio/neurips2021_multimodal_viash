#!/bin/bash

# run prior to running this script:
# bin/viash_build -q common

set -ex 

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

out_file=resources_test/common/test_resource
mkdir -p `dirname $out_file`

# remove previous output
rm $out_file*

# download raw dataset
bin/viash run src/common/datasets/download_totalvi_spleen_lymph/config.vsh.yaml -- \
  --id test_resource \
  --input "https://github.com/YosefLab/totalVI_reproducibility/raw/master/data/spleen_lymph_111.h5ad" \
  --output_rna "${out_file}.tmp.output_rna.h5ad" \
  --output_mod2 "${out_file}.tmp.output_mod2.h5ad" \
  --organism mouse

# stringent filtering to reduce the file size of test data
bin/viash run src/common/process_dataset/quality_control/config.vsh.yaml -- \
  --input_rna "${out_file}.tmp.output_rna.h5ad" \
  --input_mod2 "${out_file}.tmp.output_mod2.h5ad" \
  --output_rna "${out_file}.tmp2.output_rna.h5ad" \
  --output_mod2 "${out_file}.tmp2.output_mod2.h5ad" \
  --min_counts_per_gene 10000 \
  --min_counts_per_cell 15000 \
  --keep_genes src/common/resources/all_genes_tirosh.txt

# predetermine traintest split
bin/viash run src/common/process_dataset/split_traintest/config.vsh.yaml -- \
  --input_rna "${out_file}.tmp2.output_rna.h5ad" \
  --input_mod2 "${out_file}.tmp2.output_mod2.h5ad" \
  --output_rna "${out_file}.tmp3.output_rna.h5ad" \
  --output_mod2 "${out_file}.tmp3.output_mod2.h5ad"

# add pseudotime ordering if it's missing
bin/viash run src/common/process_dataset/pseudotime_order/config.vsh.yaml -- \
  --input_rna "${out_file}.tmp3.output_rna.h5ad" \
  --input_mod2 "${out_file}.tmp3.output_mod2.h5ad" \
  --output_rna "${out_file}.tmp4.output_rna.h5ad" \
  --output_mod2 "${out_file}.tmp4.output_mod2.h5ad"

# simulate batch if it's missing
bin/viash run src/common/process_dataset/simulate_batch/config.vsh.yaml -- \
  --input_rna "${out_file}.tmp4.output_rna.h5ad" \
  --input_mod2 "${out_file}.tmp4.output_mod2.h5ad" \
  --output_rna "${out_file}.output_rna.h5ad" \
  --output_mod2 "${out_file}.output_mod2.h5ad"
