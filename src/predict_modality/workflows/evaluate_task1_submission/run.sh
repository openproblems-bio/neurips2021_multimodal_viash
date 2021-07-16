#!/bin/bash

set +x 

NXF_VER=21.04.1 nextflow \
  run openproblems-bio/neurips2021_multimodal_viash:0.1.0 \
  -main-script src/predict_modality/workflows/evaluate_task1_method/main.nf \
  --datasets 's3://neurips2021-multimodal-public-datasets/task1_datasets/**.output_solution.h5ad'
  --predictions 'output/task1_predictions/**.h5ad' \
  --publishDir output/
  # --rootDir /home/rcannood/workspace/openproblems/neurips2021_multimodal_viash
