#!/bin/bash

echo "This script is not meant to be run! Run the commands separately, please."

TAG=1.2.0


# when building for a not-release  
bin/viash_build --max_threads 4 --config_mod '.platforms[.type == "nextflow"].separate_multiple_outputs := false'

# # build datasets for starter kits & S3
# Rscript src/common/datasets/process_inhouse_datasets/script.R
# src/joint_embedding/workflows/censor_datasets/run.sh; src/predict_modality/workflows/censor_datasets/run.sh; src/match_modality/workflows/censor_datasets/run.sh
# resources_test/run_joint_embedding.sh; resources_test/run_match_modality.sh; resources_test/run_predict_modality.sh
# src/common/create_starter_kits/create_all.sh $TAG
# aws s3 sync --delete --profile op2 output/datasets/ s3://openproblems-bio/public/phase1-data --dryrun
# aws s3 sync --delete --profile op2 output/datasets_2021-11-08/ s3://openproblems-bio/private/phase2-data/ --dryrun
# aws s3 sync --delete --profile op2 s3://openproblems-bio/private/phase2-data/phase2 s3://openproblems-bio/public/phase2-data --dryrun
# aws s3 sync --delete --profile op2 s3://openproblems-bio/private/phase2-data/phase1v2 s3://openproblems-bio/public/phase1v2-data --dryrun



rm -r target
git fetch origin
git merge origin/main

# build target folder and docker containers
bin/viash_build -m release -t $TAG --max_threads 10 --config_mod '.platforms[.type == "nextflow"].separate_multiple_outputs := false'

# run unit tests (when done right, these should all pass)
bin/viash_test -m release -t $TAG

# push docker containers to docker hub
bin/viash_push -m release -t $TAG --max_threads 10
bin/viash_push -m release -t $TAG --max_threads 4

# commit current code to release branch
git add target
git commit -m "Release $TAG"
git push 

# create new tag
git tag -a "$TAG" -m "Release $TAG"
git push --tags
