#!/bin/bash

set -ex

NXF_VER=21.04.1 nextflow \
  run openproblems-bio/neurips2021_multimodal_viash \
  -r 0.2.0 \
  -main-script src/predict_modality/workflows/run_task1_standalone_evaluation/main.nf \
  -work-dir /tmp/neurips2021_work \
  --datasets 's3://neurips2021-multimodal-public-datasets/task1_datasets/**.output_solution.h5ad' \
  --predictions 'output/task1_predictions/**.h5ad' \
  --publishDir output/