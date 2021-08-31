# Task 3: Joint Embedding

Learning of an embedded space that leverages the information of multiple modalities (e.g. for improved cell type annotation).

## Summary

The functioning of organs, tissues, and whole organisms is determined by the interplay of cells. Cells are characterised into broad types, which in turn can take on different states. Here, a cell state is made up of the sum of all processes that are occurring within the cell. We can gain insight into the state of a cell by different types of measurements: e.g., RNA expression, protein abundance, or chromatin conformation. Combining this information to describe cellular heterogeneity requires the formation of joint embeddings generated from this multimodal data. These embeddings must account for and remove possible batch effects between different measurement batches. The reward for methods that can achieve this is great: a highly resolved description of the underlying biological state of a cell that determines its function, how it interacts with other cells, and thus the cell’s role in the functioning of the whole tissue.

## Component API

### Dataset censor component

A component that censors an input datasets to the task-specific format. It expects two h5ad files containing the paired single-cell profiles using two different modalities (e.g. RNA and ADT). 

#### Input data formats

This component expects two h5ad files, `--input_mod1` and `--input_mod2`. They both contain the attributes below. If the `feature_types` of one file is `"GEX"`, then that of the other must be either `"ATAC"` or `"ADT"`.

  * `.X`: Sparse profile matrix.
  * `.obs['batch']`: Batch id.
  * `.obs['size_factors']`: The size factors computed by scran.
  * `.obs['cell_type']`: Cell type each cell belongs to.
  * `.obs['organism']`: Organism the cell was taken from (only for `"GEX"` features).
  * `.obs['S_score']`: Cell cycle score on S-phase genes (only for `"GEX"` features). Required if `var_names` are not gene symbols.
  * `.obs['G2M_score']`: Cell cycle score on G2-phase & M-phase genes (only for `"GEX"` features). Required if `var_names` are not gene symbols.
  * `.obs['pseudotime_order_GEX']`: Pseudotime values for `"GEX"` features. Required for trajectory conservation score.
  * `.obs['pseudotime_order_ATAC']`: Pseudotime values for `"ATAC"` features. Required for trajectory conservation score.
  * `.obs['pseudotime_order_ADT']`: Pseudotime values for `"ADT"` features. Required for trajectory conservation score.
  * `.var['gene_ids']`: Additional gene Ids (optional).
  * `.var['feature_types']`: Modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.
  * `.uns['organism']`: The organism of the sample. Must be one of "human", "mouse" or "synthetic".
  * `.uns['dataset_id']`: Name of the dataset.

#### Output data formats

This component should output *three* h5ad files, `--output_mod1`, `--output_mod2` and `--output_solution`. 

The `output_mod1` and `output_mod2` files contain the full profile matrices where extra metadata has been removed. These have the following attributes:

  * `.X`: Sparse count matrix.
  * `.obs['batch']`: Batch id.
  * `.obs['size_factors']`: The size factors computed by scran.
  * `.var['gene_ids']`: Additional gene Ids (optional).
  * `.var['feature_types']`: Modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.
  * `.uns['organism']`: The organism of the sample. Must be one of "human", "mouse" or "synthetic".
  * `.uns['dataset_id']`: Name of the dataset.

The `output_solution` file contains metadata on the cell profiles, which will be used to evaluate whether similar cells have been positioned closely to one another in the embedding.

  * `.X`: Sparse profile matrix of `"GEX"` features only.
  * `.obs['batch']`: Batch id.
  * `.obs['cell_type']`: Cell type each cell belongs to.
  * `.obs['size_factors']`: The size factors computed by scran.
  * `.obs['organism']`: Organism the cell was taken from (only for `"GEX"` features).
  * `.obs['S_score']`: Cell cycle score on S-phase genes (only for `"GEX"` features). Required if `var_names` are not gene symbols.
  * `.obs['G2M_score']`: Cell cycle score on G2-phase & M-phase genes (only for `"GEX"` features). Required if `var_names` are not gene symbols.
  * `.obs['pseudotime_order_GEX']`: Pseudotime values for `"GEX"` features (optional).
  * `.obs['pseudotime_order_ATAC']`: Pseudotime values for `"ATAC"` features (optional).
  * `.obs['pseudotime_order_ADT']`: Pseudotime values for `"ADT"` features (optional).
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the `"GEX"` features.
  * `.uns['organism']`: The organism of the sample. Must be one of "human", "mouse" or "synthetic".
  * `.uns['dataset_id']`: Name of the dataset.

### Method component

A component that embeds both modalities in a single embedding.

#### Input data formats

This component expects two h5ad files, `--input_mod1` and `--input_mod2`, containing the full profile matrices where extra metadata has been removed. These have the following attributes:

  * `.X`: Sparse count matrix.
  * `.obs['batch']`: Batch ids for all concatenated dataset batches.
  * `.obs['size_factors']`: The size factors computed by scran.
  * `.var['gene_ids']`: Additional gene Ids (optional).
  * `.var['feature_types']`: Modality of this file, should be equal to `"GEX"`, `"ATAC"` or `"ADT"`.
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.
  * `.uns['organism']`: The organism of the sample. Must be one of "human", "mouse" or "synthetic".
  * `.uns['dataset_id']`: Name of the dataset.

#### Output data formats

This component should output *one* h5ad file, `--output`, containing an embedding of the cells.

  * `.X`: Dense embedding matrix of the cells (at most 100 columns).
  * `.uns['dataset_id']`: Name of the dataset.
  * `.uns['method_id']`: Name of the prediction method.
  * `.obs_names`: Ids for the cells.

The embedding should have **at most 100 columns**.

### Metric component

A component which compares the predicted embedding against the ground-truth cell type information.

#### Input data formats

This component expects two h5ad files, `--input_prediction` and `--input_solution`.

The `input_prediction` file has the following attributes:

  * `.X`: Dense embedding matrix of the cells (at most 100 columns).
  * `.uns['dataset_id']`: Name of the dataset.
  * `.uns['method_id']`: Name of the prediction method.
  * `.obs_names`: Ids for the cells.

The `input_solution` file has the following attributes.

  * `.X`: Sparse profile matrix of `"GEX"` features only.
  * `.obs['batch']`: Batch id.
  * `.obs['size_factors']`: The size factors computed by scran.
  * `.obs['cell_type']`: Cell type each cell belongs to.
  * `.obs['organism']`: Organism the cell was taken from (only for `"GEX"` features).
  * `.obs['S_score']`: Cell cycle score on S-phase genes (only for `"GEX"` features). Required if `var_names` are not gene symbols.
  * `.obs['G2M_score']`: Cell cycle score on G2-phase & M-phase genes (only for `"GEX"` features). Required if `var_names` are not gene symbols.
  * `.obs['pseudotime_order_GEX']`: Pseudotime values for `"GEX"` features (optional).
  * `.obs['pseudotime_order_ATAC']`: Pseudotime values for `"ATAC"` features (optional).
  * `.obs['pseudotime_order_ADT']`: Pseudotime values for `"ADT"` features (optional).
  * `.obs_names`: Ids for the cells.
  * `.var_names`: Ids for the features.
  * `.uns['organism']`: The organism of the sample. Must be one of "human", "mouse" or "synthetic".
  * `.uns['dataset_id']`: Name of the dataset.

#### Output data formats

This component should output only *one* tsv file, `--output`, containing method and dataset metadata, and metric values which can be used to evaluate the performance of the method. It has the following columns:

  * `.uns['dataset_id']`: The name of the dataset.
  * `.uns['method_id']`: The name of the prediction method (only for `input_prediction`).
  * `.uns['metric_ids']`: The names of the outputted metrics (one or multiple).
  * `.uns['metric_values']`: The values of the outputted metrics (one or multiple, same length as `metric_ids`).

In addition, each metric component should also have a TSV file named `metrics_meta.tsv` in its directory. This TSV file should contain the columns `metric_id`, `metric_min`, `metric_max`, and `metric_higherisbetter`.
