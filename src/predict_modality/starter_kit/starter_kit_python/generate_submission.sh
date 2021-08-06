#!/bin/bash

set -e

# change these parameters if need be
MAX_MEMORY="16 GB"
MAX_TIME="10m"
MAX_CPUS="2"

# dataset location
DATASET_LOC='s3://neurips2021-multimodal-public-datasets/predict_modality/dyngen_**.output_mod[12].h5ad'

# alternatively, you could choose to download the contents to a local directory first.
# DATASET_LOC='/path/to/downloaddir/predict_modality/**.output_mod[12].h5ad'

[ ! -f config.vsh.yaml ] && echo "Couldn't find 'config.vsh.yaml!" && exit 1

# TODO: add more checks
# e.g. are nextflow and viash (>=0.5.1) are on the path?

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
# change the max time, max cpu and max memory usage to suit your needs.
viash build config.vsh.yaml -o target/nextflow -p nextflow \
  -c '.functionality.name := "method"' \
  -c '.platforms[.type == "nextflow"].publish := true' \
  -c ".platforms[.type == 'nextflow'].directive_time := '$MAX_TIME'" \
  -c ".platforms[.type == 'nextflow'].directive_memory := '$MAX_MEMORY'" \
  -c ".platforms[.type == 'nextflow'].directive_cpus := '$MAX_CPUS'"


echo ""
echo "######################################################################"
echo "##                     Fetch pipeline codebase                      ##"
echo "######################################################################"

# uncomment this to run the latest build:
# PIPELINE_VERSION=main_build

PIPELINE_VERSION=0.4.0

# pulling latest version of pipeline
export NXF_VER=21.04.1
nextflow pull openproblems-bio/neurips2021_multimodal_viash -r $PIPELINE_VERSION

echo ""
echo "######################################################################"
echo "##            Generating submission files using nextflow            ##"
echo "######################################################################"

# removing previous output
[ -f output ] && rm -r output/

nextflow \
  run openproblems-bio/neurips2021_multimodal_viash \
  -r $PIPELINE_VERSION \
  -main-script src/predict_modality/workflows/generate_submission/main.nf \
  --datasets "$DATASET_LOC" \
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
  --exclude=*.DS_Store* \
  --exclude=nextflow.config

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
