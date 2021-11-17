#!/bin/bash

set -e

# change these parameters if need be
PIPELINE_VERSION="$par_pipeline_version"

# helper functions

# get_script_dir: return the path of a bash file, following symlinks
function get_script_dir {
  SOURCE="$1"
  while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
  done
  cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd
}
function aws_sync {
  CMD="$1"
  SOURCE="$2"
  DEST="$3"
  # use aws cli if installed
  if command -v aws &> /dev/null; then
    aws s3 "$CMD" --no-sign-request "$SOURCE" "$DEST"
  # else use aws docker container instead
  else
    docker run \
      --user $(id -u):$(id -g) \
      --rm -it \
      -v $(pwd)/output:/output \
      -w / \
      amazon/aws-cli \
      s3 "$CMD" --no-sign-request "$SOURCE" "$DEST"
  fi
}

# get_latest_release: get the version number of the latest release on git
function get_latest_release {
  curl --silent "https://api.github.com/repos/$1/releases/latest" | # Get latest release from GitHub api
    grep '"tag_name":' |                                            # Get tag line
    sed -E 's/.*"([^"]+)".*/\1/'                                    # Pluck JSON value
}
LATEST_RELEASE=`get_latest_release openproblems-bio/neurips2021_multimodal_viash`

# cd to root dir of starter kit
cd `get_script_dir ${BASH_SOURCE[0]}`/..

# checking environment
scripts/0_sys_checks.sh

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
  -c '.platforms[.type == "nextflow"].publish := true'


echo ""
echo "######################################################################"
echo "##                      Sync datasets from S3                       ##"
echo "######################################################################"

# don't sync data when testing the development starter kits
if [[ $PIPELINE_VERSION != "main_build" ]]; then
  VERSION_FILE="output/datasets_phase1v2/$par_task/VERSION"

  # if the data is not found or is from a previous version starter kit,
  # sync from aws to local
  if [[ ! -f $VERSION_FILE || `cat $VERSION_FILE` != $PIPELINE_VERSION ]]; then
    mkdir -p output/datasets_phase1v2/$par_task/
    mkdir -p output/datasets_phase2/$par_task/

    # use aws cli if installed
    aws_sync sync "s3://openproblems-bio/public/phase1v2-data/$par_task/" "output/datasets_phase1v2/$par_task/"
    aws_sync sync "s3://openproblems-bio/public/phase2-data/$par_task/" "output/datasets_phase2/$par_task/"
    aws_sync cp "s3://openproblems-bio/public/phase1v2-data/meta.tsv" "output/datasets_phase1v2/meta.tsv"
    aws_sync cp "s3://openproblems-bio/public/phase2-data/meta.tsv" "output/datasets_phase2/meta.tsv"

    echo "$PIPELINE_VERSION" > $VERSION_FILE
  fi
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
  --datasets 'output/datasets_phase1v2/$par_task/**.h5ad' \
  --publishDir output/predictions/$par_task/ \
  -resume \
  -latest \
  -c scripts/nextflow.config

echo ""
echo "######################################################################"
echo "##                      Creating submission zip                     ##"
echo "######################################################################"
[ -f submission_phase1v2.zip ] && rm submission_phase1v2.zip
zip -9 -r -q submission_phase1v2.zip . \
  --exclude=*.git* \
  --exclude=*.nextflow* \
  --exclude=*work* \
  --exclude=*.DS_Store* \
  --exclude=nextflow.config \
  --exclude=output/datasets*/* \
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
echo "> evalai challenge 1111 phase $par_evalai_phase submit --file submission_phase1v2.zip --large --private"
echo ""
echo "Or this command to create a public one:"
echo "> evalai challenge 1111 phase $par_evalai_phase submit --file submission_phase1v2.zip --large --public"
echo ""
echo "Good luck!"

if [ $PIPELINE_VERSION != $LATEST_RELEASE ]; then
  echo ""
  echo "######################################################################"
  echo "##                             WARNING                              ##"
  echo "######################################################################"
  echo "A newer version of this starter kit is available! Updating to the"
  echo "latest version is strongly recommended. See README.md for more info."
fi