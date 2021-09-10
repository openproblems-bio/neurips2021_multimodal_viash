# Common dataset generators

This namespace contains common dataset generators (DataGen) for generating dataset files from raw data (from https, s3 or local storage).

## Required inputs

A DataGen component does not have any required inputs. Take a look at each component's config file for more details on its inputs.

## Required outputs

Each DataGen component should output two h5ad files, `--output_rna` and `--output_mod2`. 
These AnnData files should both have the following attributes:

  * `.X`: A sparse matrix of RNA expression counts.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.var['feature_types']`: The modality of this feature. For `output_rna`, this should be equal to `"GEX"`. For `output_mod2`, this should be equal to `"ATAC"` or `"ADT"` depending on the dataset.
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.

### Optional Attributes

For certain tasks, additional attributes are needed.
The AnnData files should contain:

  * `.obs['batch']`: A batch identifier. ¹
  * `.obs['cell_type']`: A cell type. ²
  * `.obs['organism']`: Organism the cell was taken from (only for `"GEX"` features). ³
  * `.obs['S_score']`: Cell cycle score on S-phase genes (only for `"GEX"` features). ³
  * `.obs['G2M_score']`: Cell cycle score on G2-phase & M-phase genes (only for `"GEX"` features). ³
  * `.obs['pseudotime_order_GEX']`: Pseudotime values for `"GEX"` features. ⁴
  * `.obs['pseudotime_order_ATAC']`: Pseudotime values for `"ATAC"` features. ⁴
  * `.obs['pseudotime_order_ADT']`: Pseudotime values for `"ADT"` features. ⁴

¹: Used in the 'predict modality' and 'match modality' tasks to make train/test splits.
²: Used in the 'joint embedding' task to compute various metrics.
³: Used in the 'joint embedding' task to compute the cell cycle conservation metric.
⁴: Used in the 'joint embedding' task to compute the trajectory conservation metric.

## Required resource

Each DataGen component should be accompanied by a TSV containing parameters for generating datasets from public datasets.
For example, the `download_10x_dataset` component has two required input parameters, `--id` and `--input`. The 
accompanying tsv, `input.tsv`, has two columns: `id` and `input`. `id` is the name of that dataset, and `input` is a 
https url or path name to a Cellranger h5 file.
