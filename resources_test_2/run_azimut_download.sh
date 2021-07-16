#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

echo "$REPO_ROOT"

target_dir=output
out_file=resources_test_2/common

out_dir_raw=${out_file}/azimuth_raw

mkdir -p `dirname $out_file`
mkdir -p `dirname $out_dir_raw`
    
$target_dir/download_azimuth_dataset \
  --id azimuth_ref \
  --input_count "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file" \
  --input_meta "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file&file=GSE164378%5Fsc%2Emeta%2Edata%5F3P%2Ecsv%2Egz" \
  --out_raw_count ${out_dir_raw}/GSE164378_RAW.tar \
  --out_raw_meta ${out_dir_raw}/GSE164378_sc.meta.data_3P.csv.gz \
  --output_rna ${out_file}/azimuth_ref_output_rna.h5ad \
  --output_mod2 ${out_file}/azimuth_ref_output_mod2.h5ad

