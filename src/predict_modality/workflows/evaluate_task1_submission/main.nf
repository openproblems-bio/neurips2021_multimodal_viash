nextflow.enable.dsl=2

targetDir = "${params.rootDir}/target/nextflow"

include  { calculate_task1_metrics }     from "$targetDir/predict_modality_metrics/calculate_task1_metrics/main.nf"    params(params)
include  { extract_scores }              from "$targetDir/common/extract_scores/main.nf"                               params(params)

workflow {
  main:
  // todo: update to path on s3
  def solutions = 
    Channel.fromPath("/home/rcannood/workspace/openproblems/neurips2021_multimodal_viash/output/task1_datasets/**.output_solution.h5ad")
    | map { [ it.getParent().baseName, it ] } 
    // | view { [ "Solution" ] + it }
  def predictions = 
    Channel.fromPath(params.predictions)
    | map { [ it.getParent().baseName, it ] }
    // | view { [ "Prediction" ] + it }

  predictions.join(solutions)
    | view{ [ "Combined" ] + it }
    | map { [ it[0], [ input_prediction: it[1], input_solution: it[2] ], params ] }
    | calculate_task1_metrics
    // | view{ [ "METRIC", it[0], it[1] ] }
    | map { it[1] }
    | toList()
    | map { [ "task1", it, params ] }
    | extract_scores
}

