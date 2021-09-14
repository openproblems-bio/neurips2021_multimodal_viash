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

# checking environment
scripts/0_sys_checks.sh

bin/viash test -p docker config.vsh.yaml