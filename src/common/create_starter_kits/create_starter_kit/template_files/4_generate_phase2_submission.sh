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
echo "##                      Creating submission zip                     ##"
echo "######################################################################"
[ -f submission_phase2.zip ] && rm submission_phase2.zip
zip -9 -r -q submission_phase2.zip . \
  --exclude=*.git* \
  --exclude=*.nextflow* \
  --exclude=*work* \
  --exclude=*.DS_Store* \
  --exclude=nextflow.config \
  --exclude=output/* \
  --exclude=submission*zip \
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
echo "> evalai challenge 1111 phase $par_evalai_phase2 submit --file submission_phase2.zip --large --private"
echo ""
echo "Or this command to create a public one:"
echo "> evalai challenge 1111 phase $par_evalai_phase2 submit --file submission_phase2.zip --large --public"
echo ""
echo "Good luck!"
echo ""
echo "PLEASE NOTE: the number of submissions for the phase 2 leaderboard is limited."
echo "Make sure your component is working using the phase 1 leaderboard before submitting"
echo "to phase 2."

if [ $PIPELINE_VERSION != $LATEST_RELEASE ]; then
  echo ""
  echo "######################################################################"
  echo "##                             WARNING                              ##"
  echo "######################################################################"
  echo "A newer version of this starter kit is available! Updating to the"
  echo "latest version is strongly recommended. See README.md for more info."
fi