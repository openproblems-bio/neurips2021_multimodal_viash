nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"
task = "joint_embedding"

include  { asw_batch }                  from "$targetDir/${task}_metrics/asw_batch/main.nf"                    params(params)
include  { asw_label }                  from "$targetDir/${task}_metrics/asw_label/main.nf"                    params(params)
include  { nmi }                        from "$targetDir/${task}_metrics/nmi/main.nf"                          params(params)
include  { cc_cons }                    from "$targetDir/${task}_metrics/cc_cons/main.nf"                      params(params)
include  { ti_cons }                    from "$targetDir/${task}_metrics/ti_cons/main.nf"                      params(params)
include  { graph_connectivity }         from "$targetDir/${task}_metrics/graph_connectivity/main.nf"           params(params)
include  { check_format }               from "$targetDir/${task}_metrics/check_format/main.nf"                 params(params)
include  { final_scores }               from "$targetDir/${task}_results/final_scores/main.nf"                 params(params)
include  { bind_tsv_rows }              from "$targetDir/common/bind_tsv_rows/main.nf"                         params(params)
include  { getDatasetId as get_id_predictions; getDatasetId as get_id_solutions } from "$srcDir/common/workflows/anndata_utils.nf"

params.solutions = "output/public_datasets/$task/**.output_solution.h5ad"

workflow {
  main:
  def predictions = Channel.fromPath(params.predictions) | get_id_predictions
  def solutions = Channel.fromPath(params.solutions) | get_id_solutions

  // create datasets meta
  def datasetsMeta = 
    Channel.fromPath("${params.rootDir}/results/meta_datasets.tsv")
  
  // create metrics meta
  def metricsMeta = 
    Channel.fromPath("$srcDir/$task/**/metric_meta_*.tsv")
      | toList()
      | map { [ "meta", it, params ] }
      | bind_tsv_rows
      | map{ it[1] }

  solutions.join(predictions)
    | map{ [ it[0], [ input_solution: it[1], input_prediction: it[2] ] , params ] }
    | ( asw_batch & asw_label & nmi & cc_cons & ti_cons & graph_connectivity & check_format)
    | mix
    | toList()
    | map{ [ it.collect{it[1]} ] }
    | combine(metricsMeta)
    | combine(datasetsMeta)
    | map{ [ "output", [ input: it[0], metric_meta: it[1], dataset_meta: it[2] ], params ] }
    | final_scores
}

