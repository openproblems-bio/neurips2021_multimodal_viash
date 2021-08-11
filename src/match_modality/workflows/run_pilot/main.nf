nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"
task = "match_modality"

include  { baseline_dr_nn_knn }          from "$targetDir/${task}_methods/baseline_dr_nn_knn/main.nf"          params(params)
include  { baseline_procrustes_knn }     from "$targetDir/${task}_methods/baseline_procrustes_knn/main.nf"     params(params)
include  { dummy_constant }              from "$targetDir/${task}_methods/dummy_constant/main.nf"              params(params)
include  { dummy_random }                from "$targetDir/${task}_methods/dummy_random/main.nf"                params(params)
include  { dummy_identity }              from "$targetDir/${task}_methods/dummy_identity/main.nf"              params(params)
include  { calculate_auroc }             from "$targetDir/${task}_metrics/calculate_auroc/main.nf"             params(params)
include  { extract_scores }              from "$targetDir/common/extract_scores/main.nf"                       params(params)
include  { bind_tsv_rows }               from "$targetDir/common/bind_tsv_rows/main.nf"                        params(params)
include  { getDatasetId as get_id_predictions; getDatasetId as get_id_solutions } from "$srcDir/common/workflows/anndata_utils.nf"

workflow pilot_wf {
  main:

  // get input files for methods
  def inputs = 
    Channel.fromPath("output/public_datasets/$task/**.h5ad")
      | map { [ it.getParent().baseName, it ] }
      | filter { !it[1].name.contains("output_solution") && !it[1].name.contains("output_test_sol") }
      | groupTuple
      | map { id, datas -> 
        def fileMap = datas.collectEntries { [ (it.name.split(/\./)[-2].replace("output_", "input_")), it ]}
        [ id, fileMap, params ]
      }
  
  // get solutions
  def solution = 
    Channel.fromPath("output/public_datasets/$task/**.h5ad")
      | map { [ it.getParent().baseName, it ] }
      | filter { it[1].name.contains("output_solution") || it[1].name.contains("output_test_sol") }

  // for now, code needs one of these code blocks per method.
  def b0 = inputs 
    | baseline_dr_nn_knn
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_dr_nn_knn", [ input_prediction: pred, input_solution: sol ], params ]}
  def b1 = inputs 
    | baseline_procrustes_knn
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_procrustes_knn", [ input_prediction: pred, input_solution: sol ], params ]}

  def d0 = inputs 
    | dummy_constant
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_dummy_constant", [ input_prediction: pred, input_solution: sol ], params ]}
  def d1 = inputs 
    | dummy_random
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_dummy_random", [ input_prediction: pred, input_solution: sol ], params ]}
  def d2 = solution
    | map { id, input -> [ id, input, params ] } 
    | dummy_identity
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_dummy_identity", [ input_prediction: pred, input_solution: sol ], params ]}

  def predictions = b0.mix(b1, d0, d1, d2)

  // fetch dataset ids in predictions and in solutions
  def prediction_dids = predictions | map { it[1].input_prediction } | get_id_predictions
  def solution_dids = solution | map { it[1] } | get_id_solutions

  // create solutions meta
  def solutionsMeta = solution_dids
    | map{ it[0] }
    | collectFile(name: "solutions_meta.tsv", newLine: true, seed: "dataset_id")
  
  // create metrics meta
  def metricsMeta = 
    Channel.fromPath("$srcDir/$task/**/metric_meta.tsv")
      | toList()
      | map{ [ "meta", it, params ] }
      | bind_tsv_rows
      | map{ it[1] }

  // compute metrics & combine results
  predictions
    | calculate_auroc
    | toList()
    | map{ [ it.collect{it[1]} ] }
    | combine(metricsMeta)
    | combine(solutionsMeta)
    | map{ [ task, [ input: it[0], metric_meta: it[1], dataset_meta: it[2] ], params ] }
    | extract_scores
}
