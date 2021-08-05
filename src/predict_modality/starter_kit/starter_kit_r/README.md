# Predict Modality - Starter Kit for R Users

## Getting started

This starter kit assumes you have Bash, Java and Docker installed.

Install [Viash](https://viash.io/docs/getting_started/installation/) and [NextFlow](https://www.nextflow.io/index.html#GetStarted).

Run `./generate_submission.sh` and submit your results to [eval.ai](https://eval.ai/web/challenges/challenge-page/1111/submission) by uploading the `submission.zip` file (easiest) or using the evalai-cli (recommended).

## Running your method manually
You can run the code manually on the sample dataset as follows:

```bash
$ viash run config.vsh.yaml -- \
  --input_mod1 sample_data/test_resource.mod1.h5ad \
  --input_mod2 sample_data/test_resource.mod2.h5ad \
  --output test_output.h5ad
```
    Loading dependencies
    Reading h5ad files
    Performing dimensionality reduction on the mod1 values
    Run KNN regression.
    Creating output matrix
    Creating output AnnData
    Writing predictions to file

**Tip:** You can also omit the `config.vsh.yaml` in the above command. 

**Tip #2:** Run `viash run -- --help` to view the command-line interface of your component.

## Changing the method code and/or dependencies
You can adapt `script.R` however you like. All the code between the `## VIASH START` and `## VIASH END` code blocks automatically
gets removed by viash and can be used for debugging your script.

Take a look at the `config.vsh.yaml`. It contains information on which parameters your component has, and which package dependencies
are required in order to build a Docker container for your component.

**Tip:** After making changes to the components dependencies, you will need to rebuild the docker container as follows:

```bash
$ viash run -- ---setup cachedbuild
```
    [notice] Running 'docker build -t method:dev /tmp/viashsetupdocker-method-tEX78c'

**Tip #2:** You can view the dockerfile that Viash generates from the config file using the `---dockerfile` argument:
```bash
viash run -- ---dockerfile
```
    FROM dataintuitive/randpy:r4.0_py3.8_bioc3.12

    RUN Rscript -e 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")' && \
    Rscript -e 'remotes::install_cran(c("lmds", "FNN"), repos = "https://cran.rstudio.com")'

## Competition documentation

Documentation for the competition can be found at [openproblems.bio/neurips_docs](https://openproblems.bio/neurips_docs).
The documentation for Viash is available at [viash.io/docs](https://viash.io/docs).