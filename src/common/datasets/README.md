# Common dataset generators

This namespace contains common dataset generators (DataGen) for generating dataset files from raw data (from https, s3 or local storage).

## Required inputs

A DataGen component does not have any required inputs.

## Required outputs

Each DataGen component should output two h5ad files, `--output_rna` and `--output_mod2`. 

`output_rna` is an AnnData file containing the following objects:

  * `ad.X`: A sparse matrix of RNA expression counts.
  * `ad.uns['dataset_id']`: The name of the dataset.
  * `ad.uns['modality']`: The modality of this file, should always be equal to `"RNA"`.
  * `ad.obs_names`: Ids for the cells.

`output_mod2` is an AnnData file containing the following objects:

  * `ad.X`: A sparse matrix of Antibody or ATACseq counts.
  * `ad.uns['dataset_id']`: The name of the dataset.
  * `ad.uns['modality']`: The modality of this file, should be equal to `"ATAC"` or `"Antibody"`.
  * `ad.obs_names`: Ids for the cells.

## Required resource

Each DataGen component should be accompanied by a TSV containing parameters for generating datasets from public datasets.
For example, the `download_10x_dataset` component has two required input parameters, `--id` and `--input`. The 
accompanying tsv, `input.tsv`, has two columns: `id` and `input`. `id` is the name of that dataset, and `input` is a 
https url or path name to a Cellranger h5 file.