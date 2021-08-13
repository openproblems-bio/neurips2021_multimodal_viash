#!/bin/bash

# run prior to running this script:
# bin/viash_build -q common

set -ex

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# TODO: replace evalai phase numbers!

pipeline_version=main_build

bin/viash run src/common/create_starter_kits/create_starter_kit/config.vsh.yaml -- ---setup cb

bin/viash run src/common/create_starter_kits/create_starter_kit/config.vsh.yaml -- \
  --input_dir src/predict_modality/starter_kit/starter_kit_r/ \
  --task predict_modality \
  --task_name "Predict Modality" \
  --language r \
  --language_name R \
  --block_starter 'par <- list(' \
  --evalai_phase 2276 \
  --pipeline_version $pipeline_version

bin/viash run src/common/create_starter_kits/create_starter_kit/config.vsh.yaml -- \
  --input_dir src/predict_modality/starter_kit/starter_kit_python/ \
  --task predict_modality \
  --task_name "Predict Modality" \
  --language python \
  --language_name Python \
  --block_starter 'par = dict(' \
  --evalai_phase 2276 \
  --pipeline_version $pipeline_version

bin/viash run src/common/create_starter_kits/create_starter_kit/config.vsh.yaml -- \
  --input_dir src/match_modality/starter_kit/starter_kit_r/ \
  --task match_modality \
  --task_name "Match Modality" \
  --language r \
  --language_name R \
  --block_starter 'par <- list(' \
  --evalai_phase XXXX \
  --pipeline_version $pipeline_version


bin/viash run src/common/create_starter_kits/create_starter_kit/config.vsh.yaml -- \
  --input_dir src/match_modality/starter_kit/starter_kit_python/ \
  --task match_modality \
  --task_name "Match Modality" \
  --language python \
  --language_name Python \
  --block_starter 'par = dict(' \
  --evalai_phase XXXX \
  --pipeline_version $pipeline_version

bin/viash run src/common/create_starter_kits/create_starter_kit/config.vsh.yaml -- \
  --input_dir src/joint_embedding/starter_kit/starter_kit_python/ \
  --task joint_embedding \
  --task_name "Joint Embedding" \
  --language python \
  --language_name Python \
  --block_starter 'par = dict(' \
  --evalai_phase XXXX \
  --pipeline_version $pipeline_version

bin/viash run src/common/create_starter_kits/create_starter_kit/config.vsh.yaml -- \
  --input_dir src/joint_embedding/starter_kit/starter_kit_r/ \
  --task joint_embedding \
  --task_name "Joint Embedding" \
  --language r \
  --language_name R \
  --block_starter 'par <- list(' \
  --evalai_phase XXXX \
  --pipeline_version $pipeline_version

if [ $USER == "rcannood" ]; then
  echo "Moving starter kits to Downloads dir"
  COPY_DIR="$HOME/Downloads/starter_kits"
  [ -d $COPY_DIR ] && rm -r $COPY_DIR
  cp -r output/starter_kits `dirname $COPY_DIR`

  for name in $COPY_DIR/*; do
    if [ -d $name ]; then
      out_dir="$name/output"
      echo copying to $out_dir
      mkdir -p $out_dir
      cp -r output/public_datasets $out_dir
    fi
  done
fi
