#!/bin/bash

set +x 

NXF_VER=21.04.1 nextflow \
  run . \
  -main-script /home/rcannood/workspace/openproblems/neurips2021_multimodal_viash/src/predict_modality/workflows/run_task1_standalone_evaluation/main.nf \
  --predictions 'output/task1_predictions/**.h5ad' \
  --publishDir output/ \
  --rootDir /home/rcannood/workspace/openproblems/neurips2021_multimodal_viash
