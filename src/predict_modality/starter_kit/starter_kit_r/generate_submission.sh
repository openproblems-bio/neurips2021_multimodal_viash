#!/bin/bash

set -e

[ ! -f config.vsh.yaml ] && echo "Couldn't find 'config.vsh.yaml!" && exit 1
# todo: add more checks
# e.g. are nextflow and viash (>=0.5.1) are on the path?
# todo: use s3 bucket instead of local dir

echo ""
echo "######################################################################"
echo "##              Build docker executable and container               ##"
echo "######################################################################"
viash build config.vsh.yaml -o target/docker -p docker --setup cachedbuild \
  -c '.functionality.name := "method"'


echo ""
echo "######################################################################"
echo "##                      Build nextflow module                       ##"
echo "######################################################################"
viash build config.vsh.yaml -o target/nextflow -p nextflow \
  -c '.functionality.name := "method"' \
  -c '.platforms[.type == "nextflow"].publish := true' \
  -c '.platforms[.type == "nextflow"].directive_time := "10m"' \
  -c '.platforms[.type == "nextflow"].directive_memory := "16 GB"' \
  -c '.platforms[.type == "nextflow"].directive_cpus := "4"'


echo ""
echo "######################################################################"
echo "##            Generating submission files using nextflow            ##"
echo "######################################################################"
export NXF_VER=21.04.1
# nextflow drop openproblems-bio/neurips2021_multimodal_viash
[ -f output ] && rm -r output/

# use this if you downloaded the datasets to a local folder first
# dataset_loc='/path/to/downloaddir/predict_modality/**.output_mod[12].h5ad'
dataset_loc='s3://neurips2021-multimodal-public-datasets/predict_modality/**.output_mod[12].h5ad'

nextflow \
  run openproblems-bio/neurips2021_multimodal_viash \
  -r release \
  -main-script src/predict_modality/workflows/generate_submission/main.nf \
  --datasets "$dataset_loc" \
  --publishDir output/predictions/predict_modality/ \
  -resume

echo ""
echo "######################################################################"
echo "##                      Creating submission zip                     ##"
echo "######################################################################"
[ -f submission.zip ] && rm submission.zip
zip -9 -r submission.zip . \
  --exclude=*.git* \
  --exclude=*.nextflow* \
  --exclude=*work* \
  --exclude=*.DS_Store*


# print message
echo ""
echo "######################################################################"
echo "##                        Submission summary                        ##"
echo "######################################################################"
echo "Please upload your submission at the link below:"
echo "  https://eval.ai/web/challenges/challenge-page/1111/submission"
echo ""
echo "Or use the command below create a private submission:"
echo "> evalai challenge 1111 phase 2276 submit --file submission.zip --large --private"
echo ""
echo "Or this command to create a public one:"
echo "> evalai challenge 1111 phase 2276 submit --file submission.zip --large --public"
echo ""
echo "Good luck!"