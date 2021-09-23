nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"
task = "predict_modality"

include  { mse }                         from "$targetDir/${task}_metrics/mse/main.nf"                        params(params)
include  { check_format }                from "$targetDir/${task}_metrics/check_format/main.nf"               params(params)
include  { final_scores }                from "$targetDir/${task}_results/final_scores/main.nf"               params(params)
include  { bind_tsv_rows }               from "$targetDir/common/bind_tsv_rows/main.nf"                       params(params)
include  { getDatasetId as get_id_predictions; getDatasetId as get_id_solutions } from "$srcDir/common/workflows/anndata_utils.nf"

params.solutionDir = "output/datasets/$task"

workflow {
  main:
  def predictions = Channel.fromPath(params.predictions) | get_id_predictions
  def solutions = Channel.fromPath(params.solutionDir + "/**.output_test_mod2.h5ad") | get_id_solutions

  // fetch dataset ids in predictions and in solutions
  def datasetsMeta = 
    Channel.fromPath("${params.rootDir}/results/meta_datasets.tsv")
  
  // create metrics meta
  def metricsMeta = 
    Channel.fromPath("$srcDir/$task/**/metric_meta*.tsv")
      | toList()
      | map{ [ "meta_metric", it, params ] }
      | bind_tsv_rows
      | map{ it[1] }

  solutions.join(predictions)
    | map{ [ it[0], [ input_solution: it[1], input_prediction: it[2] ] , params ] }
    | (mse & check_format)
    | mix
    | toList()
    | map{ [ it.collect{it[1]} ] }
    | combine(metricsMeta)
    | combine(datasetsMeta)
    | map{ [ "output", [ input: it[0], metric_meta: it[1], dataset_meta: it[2] ], params ] }
    | final_scores
}

