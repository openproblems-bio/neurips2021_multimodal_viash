# Task 1: Predict Modality

Predicting the profiles of one modality (e.g. protein abundance) from another (e.g. mRNA expression).

## Summary

Experimental techniques to measure multiple modalities within the same single cell are increasingly becoming available. The demand for these measurements is driven by the promise to provide a deeper insight into the state of a cell. Yet, the modalities are also intrinsically linked. We know that DNA must be accessible (ATAC data) to produce mRNA (expression data), and mRNA in turn is used as a template to produce protein (protein abundance). These processes are regulated often by the same molecules that they produce: for example, a protein may bind DNA to prevent the production of more mRNA. Understanding these regulatory processes would be transformative for synthetic biology and drug target discovery. Any method that can predict a modality from another must have accounted for these regulatory processes, but the demand for multi-modal data shows that this is not trivial.


## Interfaces

### Dataset censor component

A component that censors an input datasets to the task-specific format. It expects two h5ad files containing the paired single-cell profiles using two different modalities (e.g. RNA and ADT). 

#### Input data formats

This component component should output two h5ad files, `--input_mod1` and `--input_mod2`. They both contain the attributes below. If the `feature_types` of one file is `"GEX"`, then that of the other must be either `"ATAC"` or `"ADT"`.

  * `ad.X`: Sparse profile matrix.
  * `ad.uns['dataset_id']`: The name of the dataset.
  * `ad.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `ad.obs['batch']`: A batch identifier for the sets of cells (optional). If available, it will use the batches to derive a train/test split.
  * `ad.obs_names`: Ids for the cells.
  * `ad.var_names`: Ids for the features.

#### Output data formats

This component component should output *three* h5ad files, `--output_mod1`, `--output_mod2` and `--output_solution`. These all have the following attributes:

  * `ad.X`: Sparse profile matrix.
  * `ad.uns['dataset_id']`: The name of the dataset.
  * `ad.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * **`ad.obs['group']`: Denotes whether a cell belongs to the 'train' or the 'test' set.**
  * `ad.obs_names`: Ids for the cells.
  * `ad.var_names`: Ids for the features.

The dimensions of the three h5ad files are different;

  * `output_mod1` contains the same data as `input_mod1` (though variables and observations might be ordered differently).
  * `output_mod2` contains only the `input_mod2` data of the `"train"` cells.
  * `output_solution` contains only the `input_mod2` data of the `"test"` cells.

The `output_mod1` and `output_mod2` will be passed to the downstream regression method, whereas `output_solution` will be passed to the metric.

### Method component

A component that predicts the profile data of one modality based on the profile data of another modality. This component will receive the
profile data of both modality for a group of cells called the "train" cells, and only the profile data of one modality for the "test" cells. 
It should train a model based on the information in the train cells, and make a prediction of the unseen profile data of the second modality
for the test cells.

#### Input data formats

This component should output two h5ad files, `--input_mod1` and `--input_mod2`. They both contain the attributes below. If the `feature_types` of one file is `"GEX"`, then that of the other must be either `"ATAC"` or `"ADT"`.

  * `ad.X`: Sparse profile matrix.
  * `ad.uns['dataset_id']`: The name of the dataset.
  * `ad.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `ad.obs['group']`: Denotes whether a cell belongs to the 'train' or the 'test' set.
  * `ad.obs_names`: Ids for the cells.
  * `ad.var_names`: Ids for the features.

The dimensions of these two h5ad files are different;

  * `output_mod1` contains the modality 1 data of both the `"train"` and the `"test"` cells.
  * `output_mod2` contains only modality 2 data of the `'train'` cells.

#### Output data formats

This component should output only *one* h5ad file, `--output_prediction`, containing the predicted profile values of modality 2 for the test cells. It has the following attributes:

  * `ad.X`: Sparse profile matrix.
  * `ad.uns['dataset_id']`: The name of the dataset.
  * **`ad.uns['method_id']`: The name of the prediction method.**
  * `ad.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `ad.obs['group']`: Denotes whether a cell belongs to the 'train' or the 'test' set.
  * `ad.obs_names`: Ids for the cells.
  * `ad.var_names`: Ids for the features.



### Metric component

A component which compares the predicted profile data against the ground-truth profile data and produces one or more scores. 

#### Input data formats

This component should output two h5ad files, `--input_prediction` and `--input_solution`. The former contains the predictions of the modality 2 profile data for the test cells, whereas the latter contains the ground-truth modality 2 profile data for the test cells. Both input files should have the following interface:

  * `ad.X`: Sparse profile matrix.
  * `ad.uns['dataset_id']`: The name of the dataset.
  * `ad.uns['method_id']`: The name of the prediction method (only for `input_prediction`).
  * `ad.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `ad.obs['group']`: Denotes whether a cell belongs to the 'train' or the 'test' set.
  * `ad.obs_names`: Ids for the cells.
  * `ad.var_names`: Ids for the features.

#### Output data formats

This component should output only *one* h5ad file, `--output_prediction`, containing the predicted profile values of modality 2 for the test cells. It has the following attributes:

  * `ad.uns['dataset_id']`: The name of the dataset.
  * `ad.uns['method_id']`: The name of the prediction method (only for `input_prediction`).
  * `ad.uns['metric_ids']`: The names of the outputted metrics (one or multiple).
  * `ad.uns['metric_values']`: The values of the outputted metrics (one or multiple, same length as `metric_ids`).
  * `ad.uns['metric_moreisbetter']`: Whether or not less is better, for this metric (one or multiple, same length as `metric_ids`).