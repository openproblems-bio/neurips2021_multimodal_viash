# $par_task_name - Starter Kit for $par_language_name Users

## Getting started
This starter kit assumes you have Bash, Java and Docker installed.

Install [Viash](https://viash.io/docs/getting_started/installation/) and [NextFlow](https://www.nextflow.io/index.html#GetStarted).

Run `./generate_submission.sh` and submit your results to [eval.ai](https://eval.ai/web/challenges/challenge-page/1111/submission) by uploading the `submission.zip` file (easiest) or using the evalai-cli (recommended).

## Folder structure
In order of relevance:

    config.vsh.yaml             Metadata of your method.
    script.*                    A script containing your code.
    generate_submission.sh      A helper script for running your method on the input datasets
                                and generating a submission file to upload to eval.ai.
    sample_data/                A sample dataset for testing and debugging your code.
    nextflow.config             A config file to ensure the nextflow pipeline will be able to 
                                find your method.
    .gitignore                  A simple gitignore file.

If successful, running the `./generate_submission.sh` script will generate the following files / folders:

    submission.zip              The submission files to be uploaded to eval.ai.
    target/docker/              A standalone command-line script for running your method
                                using a Docker backend.
    target/nextflow/            A NextFlow module for using your method as part of a 
                                NextFlow pipeline.
    output/predictions/         The predictions made by your method.
    work/                       Temporary data created by the NextFlow pipeline.
    .nextflow*                  Nextflow execution output.

## Running your method manually
You can run the code manually on the sample dataset as follows:

```sh
$ viash run config.vsh.yaml -- \
  --input_mod1 sample_data/test_resource.mod1.h5ad \
  --input_mod2 sample_data/test_resource.mod2.h5ad \
  --output test_output.h5ad
```

**Tip:** You can also omit the `config.vsh.yaml` in the above command. 

**Tip #2:** Run `viash run -- --help` to view the command-line interface of your component.

## Changing the method code and/or dependencies
You can adapt the script however you like. All the code between the `## VIASH START` and `## VIASH END` code blocks automatically
gets removed by viash and can be used for debugging your script.

Take a look at the `config.vsh.yaml`. It contains information on which parameters your component has, and which package dependencies
are required in order to build a Docker container for your component.

**Tip:** After making changes to the components dependencies, you will need to rebuild the docker container as follows:

```sh
$ viash run -- ---setup cachedbuild
```
    [notice] Running 'docker build -t method:dev /tmp/viashsetupdocker-method-tEX78c'

**Tip #2:** You can view the dockerfile that Viash generates from the config file using the `---dockerfile` argument:
```sh
$ viash run -- ---dockerfile
```
    FROM dataintuitive/randpy:r4.0_py3.8_bioc3.12

    RUN Rscript -e 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")' && \
    Rscript -e 'remotes::install_cran(c("lmds", "FNN"), repos = "https://cran.rstudio.com")'

## Troubleshooting
What if running your method on the sample data works but fails when applied to the submission datasets? Given the following output:

```sh
$ ./generate_submission.sh
```
    ...
    [78/9f8fc2] NOTE: Process `method:method_process (dyngen_atac_disconnected_mod2)` terminated with an error exit status (127) -- Error is ignored
    [20/330b8a] NOTE: Process `method:method_process (dyngen_atac_bifurcating_converging_mod2)` terminated with an error exit status (127) -- Error is ignored
    Completed at: 05-Aug-2021 12:13:39
    Duration    : 1m 32s
    CPU hours   : 0.3 (100% failed)
    Succeeded   : 0
    Ignored     : 32
    Failed      : 32

You can still submit your solutions to eval.ai but you will get penalized for every failed execution. 

You can check out what went wrong by looking at the tag for the process that failed (e.g. `20/330b8a`):

```sh
$ ls -1a work/20/330b8a52b1a71489e7c53c66ef686e/
```
    .command.begin
    .command.err
    .command.log
    .command.out
    .command.run
    .command.sh
    dyngen_atac_bifurcating_converging_mod2.censor_dataset.output_mod1.h5ad
    dyngen_atac_bifurcating_converging_mod2.censor_dataset.output_mod2.h5ad
    .exitcode

You can view the `stdout` and `stderr` of this process by viewing the `.command.log` and `.command.err` files respectively. You can also try to run your component manually by editing the `VIASH START` and `VIASH END` code block in your script (see below) and running through the code manually.

```r
## VIASH START
par <- list(
    input_mod1 = "work/20/330b8a52b1a71489e7c53c66ef686e/dyngen_atac_bifurcating_converging_mod2.censor_dataset.output_mod1.h5ad",
    input_mod2 = "work/20/330b8a52b1a71489e7c53c66ef686e/dyngen_atac_bifurcating_converging_mod2.censor_dataset.output_mod2.h5ad",
    output = "debug.h5ad",
    ... other parameters
)
## VIASH END
```

## Competition documentation
Documentation for the competition can be found at [openproblems.bio/neurips_docs](https://openproblems.bio/neurips_docs).
The documentation for Viash is available at [viash.io/docs](https://viash.io/docs).