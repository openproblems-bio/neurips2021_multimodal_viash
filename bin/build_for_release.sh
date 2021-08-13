#!/bin/bash

echo "This script is not meant to be run! Run the commands separately, please."

TAG=0.5.0

rm -r target
git fetch origin
git merge origin/main

# build target folder and docker containers
bin/viash_build -m release -t $TAG --num_threads 4 \
  --config_mod '.platforms[.type == "nextflow"].separate_multiple_outputs := false' \
  --config_mod '.platforms[.type == "nextflow"].directive_memory := "10GB"' \
  --config_mod '.platforms[.type == "nextflow"].directive_time := "10 m"'

# when building for a not-release  
bin/viash_build --num_threads 4 \
  --config_mod '.platforms[.type == "nextflow"].separate_multiple_outputs := false' \
  --config_mod '.platforms[.type == "nextflow"].directive_memory := "10GB"' \
  --config_mod '.platforms[.type == "nextflow"].directive_time := "10 m"'
  

# run unit tests (when done right, these should all pass)
bin/viash_test -m release -t $TAG

# push docker containers to docker hub
bin/viash_push -m release -t $TAG

# commit current code to release branch
git add target
git commit -m "Release $TAG"
git push 

# create new tag
git tag -a "$TAG" -m "Release $TAG"
git push --tags
