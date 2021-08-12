#!/bin/bash

set -e

# change these parameters if need be
MAX_MEMORY="$par_memory"
MAX_TIME="$par_time"
MAX_CPUS="$par_cpus"
PIPELINE_VERSION="$par_pipeline_version"

[ ! -f config.vsh.yaml ] && echo "Couldn't find 'config.vsh.yaml!" && exit 1

# TODO: add more checks
# e.g. are nextflow and viash (>=0.5.1) are on the path?

echo ""
echo "######################################################################"
echo "##              Build docker executable and container               ##"
echo "######################################################################"
bin/viash build config.vsh.yaml -o target/docker -p docker --setup cachedbuild \
  -c '.functionality.name := "method"'


echo ""
echo "######################################################################"
echo "##                      Build nextflow module                       ##"
echo "######################################################################"
# change the max time, max cpu and max memory usage to suit your needs.
bin/viash build config.vsh.yaml -o target/nextflow -p nextflow \
  -c '.functionality.name := "method"' \
  -c '.platforms[.type == "nextflow"].publish := true' \
  -c ".platforms[.type == 'nextflow'].directive_time := '$MAX_TIME'" \
  -c ".platforms[.type == 'nextflow'].directive_memory := '$MAX_MEMORY'" \
  -c ".platforms[.type == 'nextflow'].directive_cpus := '$MAX_CPUS'"


echo ""
echo "######################################################################"
echo "##                      Sync datasets from S3                       ##"
echo "######################################################################"

if [ ! -d output/public_datasets/$par_task/ ]; then
  mkdir -p output/public_datasets/$par_task/

  docker run \
    --user $(id -u):$(id -g) \
    --rm -it \
    -v $(pwd)/output:/output \
    amazon/aws-cli \
    s3 sync s3://neurips2021-multimodal-public-datasets/$par_task/ /output/public_datasets/$par_task/ --no-sign-request
fi

echo ""
echo "######################################################################"
echo "##            Generating submission files using nextflow            ##"
echo "######################################################################"

export NXF_VER=21.04.1

# removing previous output
[ -d output/predictions/$par_task/ ] && rm -r output/predictions/$par_task/

bin/nextflow \
  run openproblems-bio/neurips2021_multimodal_viash \
  -r $PIPELINE_VERSION \
  -main-script src/$par_task/workflows/generate_submission/main.nf \
  --datasets 'output/public_datasets/$par_task/dyngen_**.h5ad' \
  --publishDir output/predictions/$par_task/ \
  -resume \
  -latest

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
  --exclude=nextflow.config \
  --exclude=output/public_datasets/* \
  --exclude=bin/*

# print message
echo ""
echo "######################################################################"
echo "##                        Submission summary                        ##"
echo "######################################################################"
echo "Please upload your submission at the link below:"
echo "  https://eval.ai/web/challenges/challenge-page/1111/submission"
echo ""
echo "Or use the command below create a private submission:"
echo "> evalai challenge 1111 phase $par_evalai_phase submit --file submission.zip --large --private"
echo ""
echo "Or this command to create a public one:"
echo "> evalai challenge 1111 phase $par_evalai_phase submit --file submission.zip --large --public"
echo ""
echo "Good luck!"
