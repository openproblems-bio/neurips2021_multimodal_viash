# Task 1: Predict Modality


## Dataset censor component

This component censors a dataset in the common dataset format.

### Required inputs

Each censor component should output two h5ad files, `--input_rna` and `--input_mod2`. 

`input_rna` is an AnnData file containing the following objects:

  * `ad.X`: A sparse matrix of RNA expression counts.
  * `ad.uns['dataset_id']`: The name of the dataset.
  * `ad.uns['modality']`: The modality of this file, should always be equal to `"RNA"`.
  * `ad.obs_names`: Ids for the cells.

`input_mod2` is an AnnData file containing the following objects:

  * `ad.X`: A sparse matrix of Antibody or ATACseq counts.
  * `ad.uns['dataset_id']`: The name of the dataset.
  * `ad.uns['modality']`: The modality of this file, should be equal to `"ATAC"` or `"Antibody"`.
  * `ad.obs_names`: Ids for the cells.

### Required outputs

Each censor component should output two h5ad files, `--output_mod1`, `--output_mod2` and `--output_solution`. 

`output_mod1` is an AnnData file containing expression data for modality1 for both the **train and test** cells. It contains the following objects:

  * `ad.X`: A sparse matrix of count data.
  * `ad.uns['dataset_id']`: The name of the dataset.
  * `ad.uns['modality']`: The modality of this file, should be one of `"RNA"`, `"ATAC"` and `"Antibody"`.
  * `ad.obs['experiment']`: Denotes whether a cell belongs to the 'train' or the 'test' set.
  * `ad.obs_names`: Ids for the cells.

`output_mod1` is an AnnData file containing expression data for modality2, only for the **train** cells. It contains the following objects:

  * `ad.X`: A sparse matrix of count data.
  * `ad.uns['dataset_id']`: The name of the dataset.
  * `ad.uns['modality']`: The modality of this file, should be one of `"RNA"`, `"ATAC"` and `"Antibody"`.
  * `ad.obs['experiment']`: Denotes whether a cell belongs to the 'train' or the 'test' set.
  * `ad.obs_names`: Ids for the cells.

`output_mod1` is an AnnData file containing expression data for modality2, only for the **test** cells. It contains the following objects:

  * `ad.X`: A sparse matrix of RNA expression counts.
  * `ad.uns['dataset_id']`: The name of the dataset.
  * `ad.uns['modality']`: The modality of this file, should be one of `"RNA"`, `"ATAC"` and `"Antibody"`.
  * `ad.obs_names`: Ids for the cells.