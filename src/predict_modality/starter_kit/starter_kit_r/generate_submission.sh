#!/bin/bash

set +x 

[ ! -f config.vsh.yaml ] && echo "Couldn't find 'config.vsh.yaml!" && exit 1

# building docker container and docker executable
viash build config.vsh.yaml -o target/docker -p docker --setup cachedbuild \
  -c '.functionality.name := "method"'

# building nextflow module
viash build config.vsh.yaml -o target/nextflow -p nextflow \
  -c '.functionality.name := "method"' \
  -c '.platforms[.type == "nextflow"].publish := true' \
  -c '.platforms[.type == "nextflow"].directive_time := "10m"' \
  -c '.platforms[.type == "nextflow"].directive_memory := "20 GB"'

# to do: set to repository
# nextflow drop openproblems-bio/neurips2021_multimodal_viash
# nextflow \
#   run openproblems-bio/neurips2021_multimodal_viash \
#   -main-script src/predict_modality/workflows/evaluate_task1_method/main.nf \
#   --publishDir output/task1_predictions/ \
#   -resume

# running viash pipeline
NXF_VER=21.04.1 nextflow \
  run . \
  -main-script nextflow/main.nf \
  --publishDir output/task1_predictions/ \
  -resume

[ -f submission.zip ] && rm submission.zip
zip -9 -rv submission.zip . \
  --exclude=*.git* \
  --exclude=*.nextflow* \
  --exclude=*work* \
  --exclude=*.DS_Store*

echo "Done! Please upload your submission at the link below:"
echo "  https://eval.ai/web/challenges/challenge-page/1111/submission"
echo "Alternatively, if you've set up the evalai cli, use the command below create a private submission:"
echo "  evalai challenge 1111 phase 2276 submit --file submission.zip --large --private"
echo "Or this command to create a public one:"
echo "  evalai challenge 1111 phase 2276 submit --file submission.zip --large --public"
echo "Good luck!"
# check the 'submission instructions' at the bottom of this page to set up eval cli
# https://eval.ai/web/challenges/challenge-page/1111/submission
# if you've set up the evalai-cli, you can use this command to make a private submissions
# evalai challenge 1111 phase 2276 submit --file submission.zip --large --private
