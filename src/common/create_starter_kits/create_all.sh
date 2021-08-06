#!/bin/bash

viash run src/common/create_starter_kits/create_starter_kit/config.vsh.yaml -- \
  --input_dir src/predict_modality/starter_kit/starter_kit_r/ \
  --task predict_modality \
  --task_name "Predict Modality" \
  --language r \
  --language_name R \
  --evalai_phase 2276