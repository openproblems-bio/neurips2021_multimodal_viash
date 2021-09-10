# Task 2: Match Modality

Predicting which profiles from one modality resembles a profile from another.

## Summary

While joint profiling of two modalities in the same single cell is now possible, most single-cell datasets that exist measure only a single modality. These modalities complement each other in their description of cellular state. Yet, it is challenging to analyse uni-modal datasets together when they do not share observations (cells) or a common feature space (genes, proteins, or open chromatin peaks). If we could map observations to one another across modalities, it would be possible to treat separately profiled datasets in the same manner as new multi-modal sequencing data. Mapping these modalities to one another opens up the vast amount of uni-modal single-cell datasets generated in the past years to multi-modal data analysis methods.

Unlike in task 1, where the goal was to predict _all_ values of RNA or ADT from ATAC or RNA (respectively) in each cell, the goal of this task is to identify the corresponence between single-cell profiles. Because we are only interested in matching observations, the competitors are encouraged to consider feature selection to identify the representation of the input most important for matching observations.

## Component API

### Dataset censor component

A component which partially censors a multimodal dataset. First, it will use the `.obs['is_train']` label (if available) to split the cells into train and test groups. The cell profiles are anonymised and written to file. The pairing matrices are stored as separate files.

#### Input data formats

It expects two h5ad files containing the paired single-cell profiles using two different modalities (e.g. RNA and ADT), `--input_mod1` and `--input_mod2`. They both contain the attributes below. If the `feature_types` of one file is `"GEX"`, then that of the other must be either `"ATAC"` or `"ADT"`.

  * `.X`: Sparse count matrix.
  * `.obs['batch']`: Batch id.
  * `.obs['size_factors']`: The size factors computed by scran.
  * `.obs['is_train']`: Whether or not the cell is a train or a test cell.
  * `.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['organism']`: The organism of the sample. Must be one of "human", "mouse" or "synthetic".

#### Output data formats

This component outputs *six* h5ad files, namely `--output_train_mod1`, `--output_train_mod2`, `--output_train_sol`, `--output_test_mod1`, `--output_test_mod2`, `--output_test_sol`.

The `*_mod1` and `*_mod2` h5ad files contain single-cell profiles for the two modalities for which the cell have been shuffled and anonymized. These files contain the following attributes:

  * `.X`: Sparse count matrix.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.obs['batch']`: Batch id.
  * `.obs['size_factors']`: The size factors computed by scran.
  * `.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.var_names`: Ids for the features.

The `output_train_sol` and `output_test_sol` files contain sparse matrices of which mod1 profile is paired with which mod2 profile.

  * `X`: The sparse count matrix. A value of 1 in this matrix means this modality 1 profile (row) corresponds to a modality 2 profile (column).
  * `.obs_names`: Anonymized mod1 cell ids.
  * `.var_names`: Anonymized mod2 cell ids.
  * `.uns['dataset_id']`: The name of the dataset.

### Method component

A component that predicts which profiles from one modality might match with which other profiles from the other modality. 

#### Input data formats

This component expects **five** h5ad files, `--input_train_mod1`, `--input_train_mod2`, `--input_train_sol`, `--input_test_mod1`, and `--input_test_mod2`.

The `*_mod1` and `*_mod2` h5ad files contain single-cell profiles for the two modalities for which the cell have been shuffled and anonymized. These files contain the following attributes:

  * `.X`: Sparse count matrix.
  * `.obs['batch']`: Batch id.
  * `.obs['size_factors']`: The size factors computed by scran.
  * `.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.var_names`: Ids for the features.
  * `.uns['organism']`: The organism of the sample. Must be one of "human", "mouse" or "synthetic".
  * `.uns['dataset_id']`: Name of the dataset.

The `output_train_sol` and `output_test_sol` files contain sparse matrices of which mod1 profile is paired with which mod2 profile.

  * `X`: The sparse pairing matrix. A value of 1 in this matrix means this modality 1 profile (row) corresponds to a modality 2 profile (column).
  * `.uns['dataset_id']`: The name of the dataset.

#### Output data formats

This component should output only *one* h5ad file, `--output`, containing the predicted pairings of the two input datasets.

  * `X`: The sparse pairing matrix. Dimensions N×N with at most 1000×N non-zero values, where N is the number of cells in the test set.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['method_id']`: The name of the prediction method.

If `input_test_mod1` has dimensions N×P and `input_test_mod2` has dimensions N×Q, the pairing matrix must be a sparse N×N matrix containing **at most 1000×N non-zero values**. Predictions with more than 1000×N non-zero values will be rejected by the metric component.

### Metric component

A component which compares the predicted pairing matrix against the ground-truth pairing matrix and produces one or more scores. 

#### Input data formats

This component should output two h5ad files, `--input_prediction` and `--input_solution`. Both input files should have the following interface:

  * `X`: The sparse pairing matrix
  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['method_id']`: The name of the prediction method.

If the test set contains N cells, the sparse pairing matrices should have a dimensionality of N×N with at most 1000×N non-zero values.

#### Output data formats

This component should output only *one* h5ad file, `--output`, containing metric values which can be used to evaluate the performance of the method. It has the following attributes:

  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['method_id']`: The name of the prediction method (only for `input_prediction`).
  * `.uns['metric_ids']`: The names of the outputted metrics (one or multiple).
  * `.uns['metric_values']`: The values of the outputted metrics (one or multiple, same length as `metric_ids`).

In addition, each metric component should also have a TSV file named `metrics_meta.tsv` in its directory. This TSV file should contain the columns `metric_id`, `metric_min`, `metric_max`, and `metric_higherisbetter`.
