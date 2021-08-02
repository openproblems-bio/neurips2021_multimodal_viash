# neurips2021_multimodal_viash 0.4.0

## NEW FEATURES

* Task 1: Added Babel method.
* Task 2: Added NextFlow pipeline for censoring common datasets.
* Task 3: Added NextFlow pipeline for censoring common datasets.

## MINOR CHANGES

* Renamed location of datasets produced by NextFlow pipelines.

# neurips2021_multimodal_viash 0.3.0

## NEW FEATURES

* Task 2: Added dataset censoring component.
* Task 2: Added baseline methods PCA: LMDS and UMAP.
* Task 2: Added metric: Random Forest OOB Percentage of correct predictions.
* Task 3: Added baseline method: Distance optimization.
* Task 3: Added metric: AUROC and AUPR score.

## MINOR CHANGES

* Common: Extract scores component also outputs a summary tsv.
* Common: Relevel factors in `.obs` and `.var` after QC filtering.
* Task 1: Use `.obs["batch"]` to split up cells in `"train"` and `"test"`, if batch is available.
* Task 1: Make the metric component be able to evaluate many outputs at once.

## BUG FIXES

* Common: Don't forget to set dataset ids for Azimuth and TotalVI datasets.


# neurips2021_multimodal_viash 0.2.0

## NEW FEATURES

* Common: Dataset generator for 1 Azimuth dataset
* Common: Dataset generator for 2 TotalVI Spleen Lymph datasets.
* Common: Dataset generator for 3 TotalVI 10x datasets
* Task 1: R starter Kit for task 1.


# neurips2021_multimodal_viash 0.1.0

Initial release of the OpenProblems / NeurIPS 2021 Multimodal viash evaluation pipeline.

## FEATURES

### Common components
* API descriptions for the dataset generator components.
* Dataset generator for 12 10x public datasets.
* Dataset generator for 16 dyngen simulated datasets.
* Rudimentary dataset filtering component (Needs more work).
* NextFlow pipeline for generating and filtering all the datasets.

### Task 1 components
* API descriptions for task 1 censoring, method and metric components.
* Dataset censoring component.
* Three baseline methods: Random Forests, Linear Model, and KNN Regression.
* Three metrics: Pearson correlation, Spearman correlation, and RMSE.
* Nextflow pipelines for:
  - Censoring common datasets
  - Generating a submission
  - Evaluating a submission
  - Running the baseline methods and evaluating them.

### Task 2 components
* API descriptions for task 2 censoring, method and metric components.

### Task 3 components
* API descriptions for task 3 censoring, method and metric components.
* Dataset censoring component.
* Dummy PCA method.

## KNOWN ISSUES

Unit tests are not working for now. We stopped updating them in order to implement more functionality quicker.