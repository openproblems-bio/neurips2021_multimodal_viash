#!/bin/bash

set -e

## VIASH START
par_input_dir=src/predict_modality/starter_kit/starter_kit_r/
par_output_dir=output/starter_kits
par_task=predict_modality
par_language=r
par_evalai_phase=2276
par_memory="16 GB"
par_time="10m"
par_cpus="4"
par_pipeline_version="0.4.0"
par_task_name="Predict Modality"
par_language_name=R
resources_dir=src/common/create_starter_kits/create_starter_kit/
## VIASH END

output_dir="$par_output_dir/starter_kit-$par_task-$par_language"

echo remove previous results
[[ -d $output_dir ]] && rm -r $output_dir

echo create new output dir
mkdir -p $output_dir

echo copy template files
cp $resources_dir/template_files/README.md $output_dir/
cp $resources_dir/template_files/generate_submission.sh $output_dir/
cp $resources_dir/template_files/nextflow.config $output_dir/
cp $resources_dir/template_files/LICENSE $output_dir/
cp $resources_dir/template_files/.gitignore $output_dir/

echo replace terms in templates
sed -i "s#\\\$par_task_name#$par_task_name#g" $output_dir/*
sed -i "s#\\\$par_task#$par_task#g" $output_dir/*
sed -i "s#\\\$par_language_name#$par_language_name#g" $output_dir/*
sed -i "s#\\\$par_language#$par_language#g" $output_dir/*
sed -i "s#\\\$par_evalai_phase#$par_evalai_phase#g" $output_dir/*
sed -i "s#\\\$par_memory#$par_memory#g" $output_dir/*
sed -i "s#\\\$par_time#$par_time#g" $output_dir/*
sed -i "s#\\\$par_cpus#$par_cpus#g" $output_dir/*
sed -i "s#\\\$par_pipeline_version#$par_pipeline_version#g" $output_dir/*
sed -i "s#\\\$par_block_starter#$par_block_starter#g" $output_dir/*

echo run viash dockerfile
dockerfile=$(viash run $par_input_dir/config.vsh.yaml -- ---dockerfile | sed 's#^#\t#' | sed ':a;N;$!ba;s/\n/\\n/g' | sed 's#&#\\\&#g')
sed -i "s~\\\$codeblock_dockerfile~$dockerfile~g" $output_dir/*

echo copy executables
mkdir $output_dir/bin
cp $par_bin/viash $output_dir/bin/
cp $par_bin/nextflow $output_dir/bin/

echo copy scripts
cp $par_input_dir/* $output_dir

# todo: update to multisample
echo copy sample resources
mkdir -p $output_dir/sample_data/
if [[ $par_task == "predict_modality" ]]; then
  cp $resources_dir/resources_test/$par_task/test_resource.train_mod[12].h5ad $output_dir/sample_data/
  cp $resources_dir/resources_test/$par_task/test_resource.test_mod1.h5ad $output_dir/sample_data/
elif [[ $par_task == "match_modality" ]]; then
  cp $resources_dir/resources_test/$par_task/test_resource.train_*.h5ad $output_dir/sample_data/
  cp $resources_dir/resources_test/$par_task/test_resource.test_mod[12].h5ad $output_dir/sample_data/
elif [[ $par_task == "joint_embedding" ]]; then
  cp $resources_dir/resources_test/$par_task/test_resource.mod[12].h5ad $output_dir/sample_data/
fi

echo zipping starter kit
[ -f ${output_dir}.zip ] && rm ${output_dir}.zip
cd  ${output_dir} && zip -9 -r ../$(basename $output_dir).zip *

echo starter kit is done!
