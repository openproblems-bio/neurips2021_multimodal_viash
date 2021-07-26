nextflow.enable.dsl=2

targetDir = "${params.rootDir}/target/nextflow"

include  { calculate_cor }     from "$targetDir/predict_modality_metrics/calculate_cor/main.nf"    params(params)
include  { extract_scores }              from "$targetDir/common/extract_scores/main.nf"                               params(params)

params.datasets = "s3://neurips2021-multimodal-public-datasets/task1_datasets/**.output_solution.h5ad"

workflow {
  main:
  Channel.fromList([
    [ "task1", [ input_prediction: file(params.predictions), input_solution: file(params.datasets) ], params ]
  ])
    // | view { [ "DEBUG", it[0], it[1] ] }
    | calculate_cor
    | extract_scores
}

