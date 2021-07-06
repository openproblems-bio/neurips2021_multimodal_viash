This repository contains viash components for running the
NeurIPS2021-OpenProblems benchmark pipeline for evaluating Multimodal
Data Integration methods.

## Running the whole pipeline

Download NextFlow, viash and helper components by executing:

``` sh
bin/init
```

Build all components and Docker containers (might take a while the first
time around):

``` sh
bin/viash_build
```

Use NextFlow and viash components to generate all synthetic and public
datasets:

``` sh
src/common/workflows/generate_datasets/run.sh
```

Run benchmarking pipeline for task 1:

``` sh
src/predict_modality/workflows/run_benchmark/run.sh
```

## Run an individual component

You can run an individual viash component using the `viash run` command:

``` sh
viash run src/common/datasets/download_10x_dataset/config.vsh.yaml -- \
  --id pbmc_1k_protein_v3 \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5 \
  --output output.h5ad
```

Or if you already ran `bin/viash_build`:

``` sh
target/docker/common_datasets/download_10x_dataset/download_10x_dataset \
  --id pbmc_1k_protein_v3 \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5 \
  --output output.h5ad
```

Also check out the component’s help page:

``` sh
target/docker/common_datasets/download_10x_dataset/download_10x_dataset
```

## More information

For more information on how to use `viash run`, `viash build` and
`viash test`, please take a look at the documentation at
[viash.io](https://viash.io).
