#!/bin/bash

# run prior to running this script:
# bin/viash_build -q 'common|joint_embedding'

set -ex

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

in_file=resources_test/common/openproblems_bmmc_cite_starter/openproblems_bmmc_cite_starter
out_file=resources_test/joint_embedding/openproblems_bmmc_cite_starter/openproblems_bmmc_cite_starter

# remove previous output
out_dir=`dirname $out_file`
[ -d $out_dir ] && rm -r $out_dir
mkdir -p $out_dir

viash run src/joint_embedding/datasets/censor_dataset/config.vsh.yaml -- \
  --input_mod1 ${in_file}.output_rna.h5ad \
  --input_mod2 ${in_file}.output_mod2.h5ad \
  --output_mod1 ${out_file}.mod1.h5ad \
  --output_mod2 ${out_file}.mod2.h5ad \
  --output_solution ${out_file}.solution.h5ad
  
viash run src/joint_embedding/methods/baseline_lmds/config.vsh.yaml -- \
  --input_mod1 ${out_file}.mod1.h5ad \
  --input_mod2 ${out_file}.mod2.h5ad \
  --output ${out_file}.prediction.h5ad
  
viash run src/joint_embedding/metrics/nmi/config.vsh.yaml -- \
  --input_prediction ${out_file}.prediction.h5ad \
  --input_solution ${out_file}.solution.h5ad \
  --output ${out_file}.scores.h5ad

viash run src/common/extract_scores/config.vsh.yaml -- \
  --input ${out_file}.scores.h5ad \
  --metric_meta src/joint_embedding/metrics/nmi/metric_meta_nmi.tsv \
  --output ${out_file}.scores.tsv




in_file=resources_test/common/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter
out_file=resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter

# remove previous output
out_dir=`dirname $out_file`
[ -d $out_dir ] && rm -r $out_dir
mkdir -p $out_dir

viash run src/joint_embedding/datasets/censor_dataset/config.vsh.yaml -- \
  --input_mod1 ${in_file}.output_rna.h5ad \
  --input_mod2 ${in_file}.output_mod2.h5ad \
  --output_mod1 ${out_file}.mod1.h5ad \
  --output_mod2 ${out_file}.mod2.h5ad \
  --output_solution ${out_file}.solution.h5ad
  
viash run src/joint_embedding/methods/baseline_lmds/config.vsh.yaml -- \
  --input_mod1 ${out_file}.mod1.h5ad \
  --input_mod2 ${out_file}.mod2.h5ad \
  --output ${out_file}.prediction.h5ad
  
viash run src/joint_embedding/metrics/nmi/config.vsh.yaml -- \
  --input_prediction ${out_file}.prediction.h5ad \
  --input_solution ${out_file}.solution.h5ad \
  --output ${out_file}.scores.h5ad

viash run src/common/extract_scores/config.vsh.yaml -- \
  --input ${out_file}.scores.h5ad \
  --metric_meta src/joint_embedding/metrics/nmi/metric_meta_nmi.tsv \
  --output ${out_file}.scores.tsv