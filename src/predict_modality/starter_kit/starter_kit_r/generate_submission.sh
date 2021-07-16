#!/bin/bash

set +x 

[ ! -f config.vsh.yaml ] && echo "Couldn't find 'config.vsh.yaml!" && exit 1
# todo: add more checks
# e.g. are nextflow and viash (>=0.5.1) are on the path?

echo "###########################################"
echo "## Build docker executable and container ##"
echo "###########################################"
echo ""
viash build config.vsh.yaml -o target/docker -p docker --setup cachedbuild \
  -c '.functionality.name := "method"'


echo "###########################"
echo "## Build nextflow module ##"
echo "###########################"
echo ""
viash build config.vsh.yaml -o target/nextflow -p nextflow \
  -c '.functionality.name := "method"' \
  -c '.platforms[.type == "nextflow"].publish := true' \
  -c '.platforms[.type == "nextflow"].directive_time := "10m"' \
  -c '.platforms[.type == "nextflow"].directive_memory := "20 GB"'


echo "################################################"
echo "## Generating submission files using nextflow ##"
echo "################################################"
echo ""
export NXF_VER=21.04.1
# export AWS_PROFILE=op
# nextflow drop openproblems-bio/neurips2021_multimodal_viash
nextflow \
  run openproblems-bio/neurips2021_multimodal_viash \
  -r release \
  -main-script src/predict_modality/workflows/generate_task1_submission/main.nf \
  --datasets '/home/rcannood/workspace/openproblems/neurips2021_multimodal_viash/output/task1_datasets/**.output_mod[12].h5ad' \
  --publishDir output/ \
  -resume


echo "#############################"
echo "## Creating submission zip ##"
echo "#############################"
echo ""
[ -f submission.zip ] && rm submission.zip
zip -9 -rv submission.zip . \
  --exclude=*.git* \
  --exclude=*.nextflow* \
  --exclude=*work* \
  --exclude=*.DS_Store*


# print message
echo "########################"
echo "## Submission summary ##"
echo "########################"
echo ""
echo "Done! Please upload your submission at the link below:"
echo "  https://eval.ai/web/challenges/challenge-page/1111/submission"
echo "If you set up the eval.ai CLI tool, use the command below create a private submission:"
echo "  evalai challenge 1111 phase 2276 submit --file submission.zip --large --private"
echo "Or this command to create a public one:"
echo "  evalai challenge 1111 phase 2276 submit --file submission.zip --large --public"
echo "Good luck!"