# Task 2: Match Modality

Predicting which profiles from one modality resembles a profile from another.

## Summary

While joint profiling of two modalities in the same single cell is now possible, most single-cell datasets that exist measure only a single modality. These modalities complement each other in their description of cellular state. Yet, it is challenging to analyse uni-modal datasets together when they do not share observations (cells) or a common feature space (genes, proteins, or open chromatin peaks). If we could map observations to one another across modalities, it would be possible to treat separately profiled datasets in the same manner as new multi-modal sequencing data. Mapping these modalities to one another opens up the vast amount of uni-modal single-cell datasets generated in the past years to multi-modal data analysis methods.

Unlike in task 1, where the goal was to predict _all_ values of RNA or ADT from ATAC or RNA (respectively) in each cell, the goal of this task is to identify the corresponence between single-cell profiles. Because we are only interested in matching observations, the competitors are encouraged to consider feature selection to identify the representation of the input most important for matching observations.

## Component API

### Dataset censor component

A component that censors an input datasets to the task-specific format. It expects two h5ad files containing the paired single-cell profiles using two different modalities (e.g. RNA and ADT). 

#### Input data formats

This component expects two h5ad files, `--input_mod1` and `--input_mod2`. They both contain the attributes below. If the `feature_types` of one file is `"GEX"`, then that of the other must be either `"ATAC"` or `"ADT"`.

  * `.X`: Sparse profile matrix.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.

#### Output data formats

This component should output *three* h5ad files, `--output_mod1`, `--output_mod2`, and `--output_solution`. Since this is a supervised problem, the input cells have been grouped into two groups: `'train'` and `'test'`. 

The `output_mod1` and `output_mod2` files contian the full profile matrices for the two modalities. The IDs for the cells have been removed, but for the train cells a 'pairing' is given. It has the following attributes:

  * `.X`: Sparse profile matrix.
  * **`.obs['group']`: Denotes whether a cell belongs to the 'train' or the 'test' set.**
  * `.uns['dataset_id']`: The name of the dataset.
  * `.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.var_names`: Ids for the features.
  * **`.obsm['pairing']`: The sparse pairing matrix (only in `output_mod1`).** A value of 1 in this matrix means this modality 1 profile (row) corresponds to a modality 2 profile (column). This information is only given for the train cells.

The `output_solution` file contains a sparse pairing matrix for the 'test' cells.

  * `.X`: Sparse pairing matrix for the 'test' cells.
  * `.uns['dataset_id']`: The name of the dataset.

### Method component

A component that predicts which profiles from one modality might match with which other profiles from the other modality. 

#### Input data formats

This component expects two h5ad files, `--input_mod1` and `--input_mod2`, for which the rows are shuffled and anonymised. These files have the following attributes:

  * `.X`: Sparse profile matrix.
  * `.obs['group']`: Denotes whether a cell belongs to the 'train' or the 'test' set.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.var_names`: Ids for the features.
  * `.obsm['pairing']`: The sparse pairing matrix (only in `output_mod1`). A value of 1 in this matrix means this modality 1 profile (row) corresponds to a modality 2 profile (column). This information is only given for the train cells.

#### Output data formats

This component should output only *one* h5ad file, `--output`, containing the predicted pairings of the two input datasets.

  * `.X`: Predicted sparse pairing matrix for the 'test' cells.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['method_id']`: The name of the prediction method.

If `input_mod1` has dimensions N×P and `input_mod2` has dimensions N×Q, the pairing matrix must be a sparse N×N matrix containing **at most 100×N non-zero values**. Predictions with more than 100×N non-zero values will be rejected by the metric component.

### Metric component

A component which compares the predicted pairing matrix against the ground-truth pairing matrix and produces one or more scores. 

#### Input data formats

This component should output two h5ad files, `--input_prediction` and `--input_solution`. Both input files should have the following interface:

  * `.X`: Sparse pairing matrix for the 'test' cells.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['method_id']`: The name of the prediction method (only for `input_prediction`).

If the original dataset contained N profiles, the sparse pairing matrices should have a dimensionality of N×N. The solution should have exactly N non-zero values (reflecting the correct pairing) while the prediction should have at most 100×N non-zero values. 

#### Output data formats

This component should output only *one* h5ad file, `--output`, containing metric values which can be used to evaluate the performance of the method. It has the following attributes:

  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['method_id']`: The name of the prediction method (only for `input_prediction`).
  * `.uns['metric_ids']`: The names of the outputted metrics (one or multiple).
  * `.uns['metric_values']`: The values of the outputted metrics (one or multiple, same length as `metric_ids`).

In addition, each metric component should also have a TSV file named `metrics_meta.tsv` in its directory. This TSV file should contain the columns `metric_id`, `metric_min`, `metric_max`, and `metric_higherisbetter`.