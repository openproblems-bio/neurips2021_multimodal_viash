# Common dataset generators

This namespace contains common dataset generators (DataGen) for generating dataset files from raw data (from https, s3 or local storage).

## Required inputs

A DataGen component does not have any required inputs. Take a look at each component's config file for more details on its inputs.

## Required outputs

Each DataGen component should output two h5ad files, `--output_rna` and `--output_mod2`. 
These AnnData files should both have the following attributes:

  * `ad.X`: A sparse matrix of RNA expression counts.
  * `ad.uns['dataset_id']`: The name of the dataset.
  * `ad.var['feature_types']`: The modality of this feature. For `output_rna`, this should be equal to `"GEX"`. For `output_mod2`, this should be equal to `"ATAC"` or `"ADT"` depending on the dataset.
  * `ad.obs['batch']`: A batch identifier (optional). If available, this can be used downstream to make train/test splits.
  * `ad.obs['cell_type']`: A cell type (optional). If available, this dataset can be used for task 2, otherwise not.
  * `ad.obs_names`: Ids for the cells.
  * `ad.var_names`: Ids for the features.

## Required resource

Each DataGen component should be accompanied by a TSV containing parameters for generating datasets from public datasets.
For example, the `download_10x_dataset` component has two required input parameters, `--id` and `--input`. The 
accompanying tsv, `input.tsv`, has two columns: `id` and `input`. `id` is the name of that dataset, and `input` is a 
https url or path name to a Cellranger h5 file.
