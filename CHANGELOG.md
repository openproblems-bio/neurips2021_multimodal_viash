# neurips2021_multimodal_viash 0.2.0

## NEW FEATURES

### Common components
* Dataset generator for 1 Azimuth dataset
* Dataset generator for 2 TotalVI Spleen Lymph datasets.
* Dataset generator for 3 TotalVI 10x datasets

### Task 1 components
* R starter Kit for task 1.


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