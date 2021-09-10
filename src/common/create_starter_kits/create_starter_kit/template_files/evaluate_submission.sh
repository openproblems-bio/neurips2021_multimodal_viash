#!/bin/bash

set -e

# change these parameters if need be
PIPELINE_VERSION="$par_pipeline_version"

# ViashSourceDir: return the path of a bash file, following symlinks
# usage   : ViashSourceDir ${BASH_SOURCE[0]}
# $1      : Should always be set to ${BASH_SOURCE[0]}
# returns : The absolute path of the bash file
function ViashSourceDir {
  SOURCE="$1"
  while [ -h "$SOURCE" ]; do
    DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
  done
  cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd
}

# cd to root dir of starter kit
cd `ViashSourceDir ${BASH_SOURCE[0]}`/..

[ ! -f config.vsh.yaml ] && echo "Error: Couldn't find 'config.vsh.yaml!" && exit 1


echo ""
echo "######################################################################"
echo "##                      Evaluating predictions                      ##"
echo "######################################################################"

export NXF_VER=21.04.1

nextflow run \
  openproblems-bio/neurips2021_multimodal_viash \
  -r $PIPELINE_VERSION \
  -main-script src/$par_task/workflows/evaluate_submission/main.nf \
  --solutionDir 'output/datasets/$par_task' \
  --predictions 'output/predictions/$par_task/**.h5ad' \
  --publishDir 'output/evaluation/$par_task' \
  -resume \
  -latest


# print message
echo ""
echo "######################################################################"
echo "##                        Evaluation summary                        ##"
echo "######################################################################"
echo "Evaluation results are stored at 'output/evaluation/$par_task'."