#!/bin/bash

# run prior to running this script:
# bin/viash_build -q 'common|match_modality'

set -ex

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

in_file=resources_test/common/openproblems_bmmc_cite_starter/openproblems_bmmc_cite_starter
out_file=resources_test/match_modality/openproblems_bmmc_cite_starter/openproblems_bmmc_cite_starter

# remove previous output
out_dir=`dirname $out_file`
[ -d $out_dir ] && rm -r $out_dir
mkdir -p $out_dir

bin/viash run src/match_modality/datasets/censor_dataset/config.vsh.yaml -- \
  --input_mod1 ${in_file}.output_rna.h5ad \
  --input_mod2 ${in_file}.output_mod2.h5ad \
  --output_train_mod1 ${out_file}.train_mod1.h5ad \
  --output_train_mod2 ${out_file}.train_mod2.h5ad \
  --output_train_sol ${out_file}.train_sol.h5ad \
  --output_test_mod1 ${out_file}.test_mod1.h5ad \
  --output_test_mod2 ${out_file}.test_mod2.h5ad \
  --output_test_sol ${out_file}.test_sol.h5ad
  
bin/viash run src/match_modality/methods/baseline_newwave_knnr_cbf/config.vsh.yaml -- \
  --input_train_mod1 ${out_file}.train_mod1.h5ad \
  --input_train_mod2 ${out_file}.train_mod2.h5ad \
  --input_train_sol ${out_file}.train_sol.h5ad \
  --input_test_mod1 ${out_file}.test_mod1.h5ad \
  --input_test_mod2 ${out_file}.test_mod2.h5ad \
  --output ${out_file}.prediction.h5ad

bin/viash run src/match_modality/metrics/match_probability/config.vsh.yaml -- \
  --input_prediction ${out_file}.prediction.h5ad \
  --input_solution ${out_file}.test_sol.h5ad \
  --output ${out_file}.scores.h5ad

bin/viash run src/common/extract_scores/config.vsh.yaml -- \
  --input ${out_file}.scores.h5ad \
  --metric_meta src/match_modality/metrics/match_probability/metric_meta_match_probability.tsv \
  --output ${out_file}.scores.tsv



in_file=resources_test/common/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter
out_file=resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter

# remove previous output
out_dir=`dirname $out_file`
[ -d $out_dir ] && rm -r $out_dir
mkdir -p $out_dir

bin/viash run src/match_modality/datasets/censor_dataset/config.vsh.yaml -- \
  --input_mod1 ${in_file}.output_rna.h5ad \
  --input_mod2 ${in_file}.output_mod2.h5ad \
  --output_train_mod1 ${out_file}.train_mod1.h5ad \
  --output_train_mod2 ${out_file}.train_mod2.h5ad \
  --output_train_sol ${out_file}.train_sol.h5ad \
  --output_test_mod1 ${out_file}.test_mod1.h5ad \
  --output_test_mod2 ${out_file}.test_mod2.h5ad \
  --output_test_sol ${out_file}.test_sol.h5ad
  
bin/viash run src/match_modality/methods/baseline_newwave_knnr_cbf/config.vsh.yaml -- \
  --input_train_mod1 ${out_file}.train_mod1.h5ad \
  --input_train_mod2 ${out_file}.train_mod2.h5ad \
  --input_train_sol ${out_file}.train_sol.h5ad \
  --input_test_mod1 ${out_file}.test_mod1.h5ad \
  --input_test_mod2 ${out_file}.test_mod2.h5ad \
  --output ${out_file}.prediction.h5ad

bin/viash run src/match_modality/metrics/match_probability/config.vsh.yaml -- \
  --input_prediction ${out_file}.prediction.h5ad \
  --input_solution ${out_file}.test_sol.h5ad \
  --output ${out_file}.scores.h5ad

bin/viash run src/common/extract_scores/config.vsh.yaml -- \
  --input ${out_file}.scores.h5ad \
  --metric_meta src/match_modality/metrics/match_probability/metric_meta_match_probability.tsv \
  --output ${out_file}.scores.tsv