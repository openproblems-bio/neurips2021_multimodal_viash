#!/bin/bash

set -ex

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
## VIASH END

output_dir="$par_output_dir/starter_kit-$par_task-$par_language/"

echo remove previous results
[[ -d $output_dir ]] && rm -r $output_dir

echo create new output dir
mkdir -p $output_dir

echo copy template files
cp $resources_dir/template_files/* $output_dir/

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

echo copy scripts
cp $par_input_dir/* $output_dir

echo copy sample resources
mkdir -p $output_dir/sample_data/
cp $resources_dir/resources_test/$par_task/test_resource.mod[12].h5ad $output_dir/sample_data/

echo copy .gitignore file
cp $resources_dir/template_files/.gitignore $output_dir/
cp $resources_dir/template_files/nextflow.config $output_dir/

echo starter kit is done!
