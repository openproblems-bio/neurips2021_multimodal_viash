#!/bin/bash

set -ex

# github repo & s3 datasets
# NXF_VER=21.04.1 nextflow \
#   run openproblems-bio/neurips2021_multimodal_viash \
#   -r main_build \
#   -main-script src/joint_embedding/workflows/evaluate_submission/main.nf \
#   -work-dir /tmp/neurips2021_work \
#   --solutions 's3://neurips2021-multimodal-public-datasets/joint_embedding/**.output_solution.h5ad' \
#   --predictions 'output/predictions/joint_embedding/**.h5ad' \
#   --publishDir output/predictions/joint_embedding/

# local repo & local datasets
# NXF_VER=21.04.1 nextflow run \
#   . \
#   -main-script src/joint_embedding/workflows/evaluate_submission/main.nf \
#   -work-dir /tmp/neurips2021_work \
#   --solutions '/home/rcannood/workspace/openproblems/neurips2021_multimodal_viash/output/public_datasets/joint_embedding/dyngen_*/**.output_solution.h5ad' \
#   --predictions '/home/rcannood/Downloads/submission/output/predictions/joint_embedding/dyngen_*/**.h5ad' \
#   --publishDir /home/rcannood/Downloads/submission/output/predictions/joint_embedding/ \
#   -resume

# local repo & local datasets
NXF_VER=21.04.1 nextflow run \
  . \
  -main-script src/joint_embedding/workflows/evaluate_submission/main.nf \
  -work-dir /tmp/neurips2021_work \
  --solutions 'output/public_datasets/joint_embedding/**.output_solution.h5ad' \
  --predictions '/home/rcannood/Downloads/starter_kits/starter_kit-joint_embedding-r/output/predictions/joint_embedding/**.h5ad' \
  --publishDir /home/rcannood/Downloads/starter_kits/starter_kit-joint_embedding-r/output \
  -resume