#!/bin/bash

echo "This script is not meant to be run! Run the commands separately, please."

TAG=0.2.0

# build target folder and docker containers
bin/viash_build -m release -t $TAG

# run unit tests (when done right, these should all pass)
bin/viash_test -m release -t $TAG

# push docker containers to docker hub
bin/viash_push -m release -t $TAG

# commit current code to release branch
git commit -m "Release $TAG"
git push 

# create new tag
git tag -a "Release $TAG"
git push --tags

