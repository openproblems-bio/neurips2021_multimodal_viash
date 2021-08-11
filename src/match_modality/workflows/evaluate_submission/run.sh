#!/bin/bash

set -ex

# github repo & s3 datasets
NXF_VER=21.04.1 nextflow \
  run openproblems-bio/neurips2021_multimodal_viash \
  -r main_build \
  -main-script src/match_modality/workflows/evaluate_submission/main.nf \
  -work-dir /tmp/neurips2021_work \
  --solutions 's3://neurips2021-multimodal-public-datasets/match_modality/**.output_test_sol.h5ad' \
  --predictions 'output/predictions/match_modality/**.h5ad' \
  --publishDir output/predictions/match_modality/

# local repo & local datasets
# NXF_VER=21.04.1 nextflow run \
#   . \
#   -main-script src/match_modality/workflows/evaluate_submission/main.nf \
#   -work-dir /tmp/neurips2021_work \
#   --solutions '/home/rcannood/workspace/openproblems/neurips2021_multimodal_viash/output/public_datasets/match_modality/dyngen_*/**.output_test_sol.h5ad' \
#   --predictions '/home/rcannood/Downloads/submission/output/predictions/match_modality/dyngen_*/**.h5ad' \
#   --publishDir /home/rcannood/Downloads/submission/output/predictions/match_modality/ \
#   -resume