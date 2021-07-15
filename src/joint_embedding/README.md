# Task 2: Joint Embedding

Learning of an embedded space that leverages the information of multiple modalities (e.g. for improved cell type annotation).

## Summary

The functioning of organs, tissues, and whole organisms is determined by the interplay of cells. Cells are characterised into broad types, which in turn can take on different states. Here, a cell state is made up of the sum of all processes that are occurring within the cell. We can gain insight into the state of a cell by different types of measurements: e.g., RNA expression, protein abundance, or chromatin conformation. Combining this information to describe cellular heterogeneity requires the formation of joint embeddings generated from this multimodal data. These embeddings must account for and remove possible batch effects between different measurement batches. The reward for methods that can achieve this is great: a highly resolved description of the underlying biological state of a cell that determines its function, how it interacts with other cells, and thus the cell’s role in the functioning of the whole tissue.

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

This component should output *three* h5ad files, `--output_mod1`, `--output_mod2` and `--output_solution`. 

The `output_mod1` and `output_mod2` files contain the full profile matrices where extra metadata has been removed. These have the following attributes:

  * `.X`: Sparse profile matrix.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.

The `output_solution` file contains metadata on the cell profiles, which will be used to evaluate whether similar cells have been positioned closely to one another in the embedding.

  * `.obs["cell_type"]`: The cell type each cell belongs to.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.obs_names`: Ids for the cells.

### Method component

A component that embeds both modalities in a single embedding.

#### Input data formats

This component expects two h5ad files, `--input_mod1` and `--input_mod2`, containing the full profile matrices where extra metadata has been removed. These have the following attributes:

  * `.X`: Sparse profile matrix.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.var['feature_types']`: The modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.

#### Output data formats

This component should output *one* h5ad file, `--output_prediction`, containing an embedding of the cells.

  * `.X`: The embedding matrix of the cells.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['method_id']`: The name of the prediction method.
  * `.obs_names`: Ids for the cells.

The embedding should have **at most 10 columns**.

### Metric component

A component which compares the predicted embedding against the ground-truth cell type information.

#### Input data formats

This component should output two h5ad files, `--input_prediction` and `--input_solution`.

The `input_prediction` file has the following attributes:

  * `.X`: The embedding matrix of the cells (at most 10 columns).
  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['method_id']`: The name of the prediction method.
  * `.obs_names`: Ids for the cells.

The `input_solution` file has the following attributes.

  * `.obs["cell_type"]`: The cell type each cell belongs to.
  * `.uns['dataset_id']`: The name of the dataset.
  * `.obs_names`: Ids for the cells.

#### Output data formats

This component should output only *one* h5ad file, `--output`, containing metric values which can be used to evaluate the performance of the method. It has the following attributes:

  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['method_id']`: The name of the prediction method (only for `input_prediction`).
  * `.uns['metric_ids']`: The names of the outputted metrics (one or multiple).
  * `.uns['metric_values']`: The values of the outputted metrics (one or multiple, same length as `metric_ids`).
  * `.uns['metric_moreisbetter']`: Whether or not less is better, for this metric (one or multiple, same length as `metric_ids`).