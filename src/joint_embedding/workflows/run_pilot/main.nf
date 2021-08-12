nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"
task = "joint_embedding"

include  { baseline_lmds }              from "$targetDir/${task}_methods/baseline_lmds/main.nf"                params(params)
include  { baseline_pca }               from "$targetDir/${task}_methods/baseline_pca/main.nf"                 params(params)
include  { baseline_umap }              from "$targetDir/${task}_methods/baseline_umap/main.nf"                params(params)
include  { dummy_random }               from "$targetDir/${task}_methods/dummy_random/main.nf"                 params(params)
include  { dummy_zeros }                from "$targetDir/${task}_methods/dummy_zeros/main.nf"                  params(params)
include  { totalvi }                    from "$targetDir/${task}_methods/totalvi/main.nf"                      params(params)
include  { calculate_rf_oob }           from "$targetDir/${task}_metrics/calculate_rf_oob/main.nf"             params(params)
include  { calculate_totalvi_metrics }  from "$targetDir/${task}_metrics/calculate_totalvi_metrics/main.nf"    params(params)
include  { ari }                        from "$targetDir/${task}_metrics/ari/main.nf"                          params(params)
include  { asw_batch }                  from "$targetDir/${task}_metrics/asw_batch/main.nf"                    params(params)
include  { asw_label }                  from "$targetDir/${task}_metrics/asw_label/main.nf"                    params(params)
include  { nmi }                        from "$targetDir/${task}_metrics/nmi/main.nf"                          params(params)
include  { check_format }               from "$targetDir/${task}_metrics/check_format/main.nf"                 params(params)
include  { extract_scores }             from "$targetDir/common/extract_scores/main.nf"                        params(params)
include  { bind_tsv_rows }              from "$targetDir/common/bind_tsv_rows/main.nf"                         params(params)
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
    | baseline_lmds
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_lmds", [ input_prediction: pred, input_solution: sol ], params ]}
  def b1 = inputs 
    | baseline_pca
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_pca", [ input_prediction: pred, input_solution: sol ], params ]}
  def b2 = inputs 
    | baseline_umap
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_umap", [ input_prediction: pred, input_solution: sol ], params ]}

  def d0 = inputs 
    | dummy_random
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_dummy_random", [ input_prediction: pred, input_solution: sol ], params ]}
  def d1 = inputs 
    | dummy_zeros
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_dummy_zeros", [ input_prediction: pred, input_solution: sol ], params ]}

  def m0 = inputs 
    | totalvi
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_totalvi", [ input_prediction: pred, input_solution: sol ], params ]}

  def predictions = b0.mix(b1, b2, d0, d1, m0)

  // fetch dataset ids in predictions and in solutions
  def prediction_dids = predictions | map { it[1].input_prediction } | get_id_predictions
  def solution_dids = solution | map { it[1] } | get_id_solutions

  // create solutions meta
  def solutionsMeta = solution_dids
    | map{ it[0] }
    | collectFile(name: "solutions_meta.tsv", newLine: true, seed: "dataset_id")
  
  // create metrics meta
  def metricsMeta = 
    Channel.fromPath("$srcDir/$task/**/metric_meta_*.tsv")
      | toList()
      | map { [ "meta", it, params ] }
      | bind_tsv_rows
      | map{ it[1] }

  // compute metrics & combine results
  predictions
    | (calculate_rf_oob & calculate_totalvi_metrics & ari & asw_batch & asw_label & nmi)
    | mix
    | toList()
    | map{ [ it.collect{it[1]} ] }
    | combine(metricsMeta)
    | combine(solutionsMeta)
    | map{ [ "output", [ input: it[0], metric_meta: it[1], dataset_meta: it[2] ], params ] }
    | extract_scores
}
