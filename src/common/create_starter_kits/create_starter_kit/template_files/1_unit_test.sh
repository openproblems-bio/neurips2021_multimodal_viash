#!/bin/bash

set -e

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
# cd to root dir of starter kit
cd `get_script_dir ${BASH_SOURCE[0]}`/..

[ ! -f config.vsh.yaml ] && echo "Error: Couldn't find 'config.vsh.yaml!" && exit 1


# check docker availability
if ! docker info > /dev/null 2>&1; then
  echo "Docker doesn't seem to be running. Try running 'docker run hello-world'."
  exit 1
fi

bin/viash test -p docker config.vsh.yaml