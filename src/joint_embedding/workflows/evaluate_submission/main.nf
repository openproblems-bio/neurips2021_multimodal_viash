nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"

include  { calculate_rf_oob }           from "$targetDir/joint_embedding_metrics/calculate_rf_oob/main.nf"             params(params)
include  { calculate_totalVI_metrics }  from "$targetDir/joint_embedding_metrics/calculate_totalVI_metrics/main.nf"    params(params)
include  { ari }                        from "$targetDir/joint_embedding_metrics/ari/main.nf"                          params(params)
include  { asw_batch }                  from "$targetDir/joint_embedding_metrics/asw_batch/main.nf"                    params(params)
include  { asw_label }                  from "$targetDir/joint_embedding_metrics/asw_label/main.nf"                    params(params)
include  { nmi }                        from "$targetDir/joint_embedding_metrics/nmi/main.nf"                          params(params)
include  { extract_scores }             from "$targetDir/common/extract_scores/main.nf"                                params(params)
include  { bind_tsv_rows }              from "$targetDir/common/bind_tsv_rows/main.nf"                                 params(params)
include  { getDatasetId as get_id_predictions; getDatasetId as get_id_solutions } from "$srcDir/common/workflows/anndata_utils.nf"

params.solutions = "s3://neurips2021-multimodal-public-datasets/joint_embedding/**.output_solution.h5ad"

workflow {
  main:
  def predictions = Channel.fromPath(params.predictions) | get_id_predictions
  def solutions = Channel.fromPath(params.solutions) | get_id_solutions

  // create solutions meta
  def solutionsMeta = solutions
    | map{ it[0] }
    | collectFile(name: "solutions_meta.tsv", newLine: true, seed: "dataset_id")
  
  // create metrics meta
  def metricsMeta = 
    Channel.fromPath("$srcDir/joint_embedding/**/metric_meta.tsv")
      | toList()
      | map{ [ "meta", it, params ] }
      | bind_tsv_rows
      | map{ it[1] }

  solutions.join(predictions)
    | map{ [ it[0], [ input_solution: it[1], input_prediction: it[2] ] , params ] }
    | (calculate_rf_oob & calculate_totalVI_metrics & ari & asw_batch & asw_label & nmi)
    | mix
    | toList()
    | map{ [ it.collect{it[1]} ] }
    | combine(metricsMeta)
    | combine(solutionsMeta)
    | map{ [ "output", [ input: it[0], metric_meta: it[1], dataset_meta: it[2] ], params ] }
    | extract_scores
}
