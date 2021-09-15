#!/bin/bash

echo "This script is not meant to be run! Run the commands separately, please."

TAG=1.0.2

rm -r target
git fetch origin
git merge origin/main

# build target folder and docker containers
bin/viash_build -m release -t $TAG --max_threads 4 --config_mod '.platforms[.type == "nextflow"].separate_multiple_outputs := false'

# when building for a not-release  
bin/viash_build --max_threads 4 --config_mod '.platforms[.type == "nextflow"].separate_multiple_outputs := false'
  

# run unit tests (when done right, these should all pass)
bin/viash_test -m release -t $TAG

# push docker containers to docker hub
bin/viash_push -m release -t $TAG --max_threads 10
bin/viash_push -m release -t $TAG --max_threads 4

# build starter kits
src/common/create_starter_kits/create_all.sh $TAG

# commit current code to release branch
git add target
git commit -m "Release $TAG"
git push 

# create new tag
git tag -a "$TAG" -m "Release $TAG"
git push --tags

aws s3 sync --delete --profile op2 output/datasets/ s3://openproblems-bio/public/phase1-data --dryrun
