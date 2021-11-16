#!/bin/bash

set -e

## VIASH START
par_src=src/
par_output_dir=output/starter_kits
par_task=predict_modality
par_language=r
par_evalai_phase=2276
par_evalai_phase2=2276
par_memory="16 GB"
par_time="10m"
par_cpus="4"
par_pipeline_version="0.4.0"
par_task_name="Predict Modality"
par_language_name=R
resources_dir=src/common/create_starter_kits/create_starter_kit/
## VIASH END

input_dir="$par_src/$par_task/starter_kit/starter_kit_$par_language"
output_dir="$par_output_dir/starter_kit-$par_task-$par_language"

echo ">> Creating $par_task_name starter kit for $par_language_name users."

echo "  Remove previous results"
[[ -d $output_dir ]] && rm -r $output_dir

echo "  Create new output dir"
mkdir -p $output_dir/scripts

echo "  Copy template files"
cp $resources_dir/template_files/README.md $output_dir/
cp $resources_dir/template_files/0_sys_checks.sh $output_dir/scripts/
cp $resources_dir/template_files/1_unit_test.sh $output_dir/scripts/
cp $resources_dir/template_files/2_generate_submission.sh $output_dir/scripts/
cp $resources_dir/template_files/3_evaluate_submission.sh $output_dir/scripts/
cp $resources_dir/template_files/4_generate_phase2_submission.sh $output_dir/scripts/
cp $resources_dir/template_files/nextflow.config $output_dir/scripts/
cp $resources_dir/template_files/LICENSE $output_dir/
cp $resources_dir/template_files/.gitignore $output_dir/

echo "  Run viash dockerfile"
dockerfile=$(viash run $input_dir/config.vsh.yaml -- ---dockerfile | sed 's#^#\t#' | sed ':a;N;$!ba;s/\n/\\n/g' | sed 's#&#\\\&#g')

echo "  Replace terms in templates"
for file in $(find $output_dir/ -type f); do
  sed -i "s#\\\$par_task_name#$par_task_name#g" $file
  sed -i "s#\\\$par_task#$par_task#g" $file
  sed -i "s#\\\$par_language_name#$par_language_name#g" $file
  sed -i "s#\\\$par_language#$par_language#g" $file
  sed -i "s#\\\$par_evalai_phase2#$par_evalai_phase2#g" $file
  sed -i "s#\\\$par_evalai_phase#$par_evalai_phase#g" $file
  sed -i "s#\\\$par_memory#$par_memory#g" $file
  sed -i "s#\\\$par_time#$par_time#g" $file
  sed -i "s#\\\$par_cpus#$par_cpus#g" $file
  sed -i "s#\\\$par_pipeline_version#$par_pipeline_version#g" $file
  sed -i "s~\\\$codeblock_dockerfile~$dockerfile~g" $file
done


echo "  Copy executables"
mkdir $output_dir/bin
cp $par_bin/viash $output_dir/bin/
cp $par_bin/nextflow $output_dir/bin/

echo "  Copy scripts"
cp $input_dir/config.vsh.yaml $output_dir/config.vsh.yaml
cat $input_dir/script.$par_language_ext  | sed "s#resources_test/$par_task#sample_data#g" > $output_dir/script.$par_language_ext

unit_test=$par_src/$par_task/unit_tests/test_method.$par_language_ext
if [[ -f $unit_test ]]; then
  echo "  Copy unit test"
  cat $unit_test | sed "s#resources_test/$par_task#sample_data#g" > $output_dir/test.$par_language_ext
fi

echo "  Build check_format component"
viash build $par_src/$par_task/metrics/check_format/config.vsh.yaml -p docker -o $output_dir/bin/

# todo: update to multisample
echo "  Copy sample resources"
rsync -avzr $resources_dir/resources_test/$par_task/ $output_dir/sample_data/ \
  --include="*/" --include="*mod[12].h5ad" --include="*solution.h5ad" --include="*sol.h5ad" --exclude="*" # --dry-run

echo "  Zipping starter kit"
[ -f ${output_dir}.zip ] && rm ${output_dir}.zip
cd  ${output_dir} && zip -9 -q -r ../$(basename $output_dir).zip *

echo "  Starter kit is done!"
