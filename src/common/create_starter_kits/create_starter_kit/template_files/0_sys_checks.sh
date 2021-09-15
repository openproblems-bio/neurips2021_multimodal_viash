#!/bin/bash

set -e

# change these parameters if need be
PIPELINE_VERSION="$par_pipeline_version"

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

# cd to root dir of starter kit
cd `get_script_dir ${BASH_SOURCE[0]}`/..

if [ ! -f config.vsh.yaml ]; then
  echo "Error: Couldn't find 'config.vsh.yaml!"
  exit 1
fi

LATEST_RELEASE=`get_latest_release openproblems-bio/neurips2021_multimodal_viash`
if [ $PIPELINE_VERSION == "main_build" ]; then
  echo "Warning: Running dev build pipelines. Only do this if you know what you're doing!"
elif [ $PIPELINE_VERSION != $LATEST_RELEASE ]; then
  echo "Warning: A newer version of this starter kit is available! Updating to"
  echo "the latest version is strongly recommended. See README.md for more info."
  echo "Continuing in 10 seconds."
  sleep 10
fi

# check docker availability
if ! docker info > /dev/null 2>&1; then
  echo "Error: Docker does not seem to be running. Try running 'docker run hello-world'."
  exit 1
fi

if type -p java > /dev/null; then
  :
elif [[ -n "$JAVA_HOME" ]] && [[ -x "$JAVA_HOME/bin/java" ]]; then
  :
else
  echo "Error: Java was not found."
  exit 1
fi