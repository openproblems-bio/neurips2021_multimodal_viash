nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"

include  { calculate_cor }             from "$targetDir/predict_modality_metrics/calculate_cor/main.nf"    params(params)
include  { extract_scores }            from "$targetDir/common/extract_scores/main.nf"                     params(params)
include  { bind_tsv_rows }             from "$targetDir/common/bind_tsv_rows/main.nf"                      params(params)
include  { getDatasetId as DID0; getDatasetId as DID1 } from "$srcDir/common/workflows/anndata_utils.nf"

params.solutions = "s3://neurips2021-multimodal-public-datasets/predict_modality/**.output_solution.h5ad"

workflow {
  main:
  def predictions = Channel.fromPath(params.predictions) | DID0
  def solutions = Channel.fromPath(params.solutions) | DID1

  // create solutions meta
  def solutionsMeta = solutions
    | map{ it[0] }
    | collectFile(name: "solutions_meta.tsv", newLine: true, seed: "dataset_id")
  
  // create metrics meta
  def metricsMeta = 
    Channel.fromPath("$srcDir/predict_modality/**/metric_meta.tsv")
      | toList()
      | map{ [ "meta", it, params ] }
      | bind_tsv_rows
      | map{ it[1] }

  solutions.join(predictions)
    | map{ [ it[0], [ input_solution: it[1], input_prediction: it[2] ] , params ] }
    | calculate_cor
    | toList()
    | map{ [ it.collect{it[1]} ] }
    | combine(metricsMeta)
    | combine(solutionsMeta)
    | map{ [ "output", [ input: it[0], metric_meta: it[1], dataset_meta: it[2] ], params ] }
    | extract_scores
}

