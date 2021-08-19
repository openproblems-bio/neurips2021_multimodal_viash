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

bin/viash run src/common/datasets/simulate_dyngen_dataset/config.vsh.yaml -- \
  --id test_resource \
  --num_cells 200 \
  --num_genes 150 \
  --num_simulations 20 \
  --store_protein \
  --num_threads 10 \
  --seed 1 \
  --output_rna "${out_file}.tmp.output_rna.h5ad" \
  --output_mod2 "${out_file}.tmp.output_mod2.h5ad" \
  --plot "${out_file}.plot.pdf"

# stringent filtering to reduce the file size of test data
bin/viash run src/common/process_dataset/quality_control/config.vsh.yaml -- \
  --input_rna "${out_file}.tmp.output_rna.h5ad" \
  --input_mod2 "${out_file}.tmp.output_mod2.h5ad" \
  --output_rna "${out_file}.tmp2.output_rna.h5ad" \
  --output_mod2 "${out_file}.tmp2.output_mod2.h5ad" \
  --min_counts_per_gene 10000 \
  --min_counts_per_cell 15000

# predetermine traintest split
bin/viash run src/common/process_dataset/split_traintest/config.vsh.yaml -- \
  --input_rna "${out_file}.tmp2.output_rna.h5ad" \
  --input_mod2 "${out_file}.tmp2.output_mod2.h5ad" \
  --output_rna "${out_file}.output_rna.h5ad" \
  --output_mod2 "${out_file}.output_mod2.h5ad"

rm ${out_file}.tmp*
