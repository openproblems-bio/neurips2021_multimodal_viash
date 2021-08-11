# Task 1: Predict Modality

Predicting the profiles of one modality (e.g. protein abundance) from another (e.g. mRNA expression).

## Summary

Experimental techniques to measure multiple modalities within the same single cell are increasingly becoming available. The demand for these measurements is driven by the promise to provide a deeper insight into the state of a cell. Yet, the modalities are also intrinsically linked. We know that DNA must be accessible (ATAC data) to produce mRNA (expression data), and mRNA in turn is used as a template to produce protein (protein abundance). These processes are regulated often by the same molecules that they produce: for example, a protein may bind DNA to prevent the production of more mRNA. Understanding these regulatory processes would be transformative for synthetic biology and drug target discovery. Any method that can predict a modality from another must have accounted for these regulatory processes, but the demand for multi-modal data shows that this is not trivial.


## Component API

### Dataset censor component

A component which partially censors a multimodal dataset, using the `.obs['is_train']` label (if available) to split the cells into train and test groups.

#### Input data formats

This component expects two h5ad files, `--input_mod1` and `--input_mod2`. They both contain the attributes below. If the `feature_types` of one file is `"GEX"`, then that of the other must be either `"ATAC"` or `"ADT"`.

  * `.X`: Sparse profile matrix.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.obs['batch']`: A batch identifier for the sets of cells (optional). If available, it will use the batches to derive a train/test split.
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.

#### Output data formats

This component outputs *four* h5ad files, namely `--output_train_mod1`, `--output_train_mod2`, `--output_test_mod1`, and `--output_test_mod2`. Each of these files contain the following attributes:

  * `.X`: Sparse profile matrix.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.

The dimensions of the four h5ad files are different;

  * `output_train_mod1` contains the cells in `input_mod1` for which `.obs['is_train']` is true.
  * `output_train_mod2` contains the cells in `input_mod2` for which `.obs['is_train']` is true.
  * `output_test_mod1` contains the cells in `input_mod1` for which `.obs['is_train']` is false.
  * `output_test_mod2` contains the cells in `input_mod2` for which `.obs['is_train']` is false.

The train and test mod1 files will be passed to the downstream regression method, whereas `output_test_mod2` will be passed to the metric.

### Method component

A component that predicts the profile data of one modality based on the profile data of another modality. This component will receive the
profile data of both modality for a group of cells called the "train" cells, and only the profile data of one modality for the "test" cells. 
It should train a model based on the information in the train cells, and make a prediction of the unseen profile data of the second modality
for the test cells.

#### Input data formats

This component expects three inputs, `--input_train_mod1`, `--input_train_mod2`, and `--input_test_mod1`. These should contain the attributes below. 

  * `.X`: Sparse profile matrix.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.

#### Output data formats

This component should output only *one* h5ad file, `--output`, containing the predicted profile values of modality 2 for the test cells. It has the following attributes:

  * `.X`: Sparse profile matrix.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['method_id']`: The name of the prediction method.
  * `.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.


### Metric component

A component which compares the predicted profile data against the ground-truth profile data and produces one or more scores. 

#### Input data formats

This component expects two h5ad files, `--input_prediction` and `--input_solution`. The former contains the predictions of the modality 2 profile data for the test cells, whereas the latter contains the ground-truth modality 2 profile data for the test cells. Both input files should have the following interface:

  * `.X`: Sparse profile matrix.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['method_id']`: The name of the prediction method (only for `input_prediction`).
  * `.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.

#### Output data formats

This component should output only *one* h5ad file, `--output`, containing metric values which can be used to evaluate the performance of the method. It has the following attributes:

  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['method_id']`: The name of the prediction method (only for `input_prediction`).
  * `.uns['metric_ids']`: The names of the outputted metrics (one or multiple).
  * `.uns['metric_values']`: The values of the outputted metrics (one or multiple, same length as `metric_ids`).

In addition, each metric component should also have a TSV file named `metrics_meta.tsv` in its directory. This TSV file should contain the columns `metric_id`, `metric_min`, `metric_max`, and `metric_higherisbetter`.
