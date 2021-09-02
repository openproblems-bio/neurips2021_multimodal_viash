nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"
task = "joint_embedding"

include  { baseline_lmds }              from "$targetDir/${task}_methods/baseline_lmds/main.nf"                params(params)
include  { baseline_pca }               from "$targetDir/${task}_methods/baseline_pca/main.nf"                 params(params)
include  { baseline_mnn }               from "$targetDir/${task}_methods/baseline_mnn/main.nf"                 params(params)
include  { baseline_umap }              from "$targetDir/${task}_methods/baseline_umap/main.nf"                params(params)
include  { baseline_totalvi }           from "$targetDir/${task}_methods/baseline_totalvi/main.nf"             params(params)
include  { baseline_newwave }           from "$targetDir/${task}_methods/baseline_newwave/main.nf"             params(params)
include  { dummy_random }               from "$targetDir/${task}_methods/dummy_random/main.nf"                 params(params)
include  { dummy_zeros }                from "$targetDir/${task}_methods/dummy_zeros/main.nf"                  params(params)
include  { dummy_solution }             from "$targetDir/${task}_methods/dummy_solution/main.nf"               params(params)
include  { rfoob }                      from "$targetDir/${task}_metrics/rfoob/main.nf"                        params(params)
include  { latent_mixing }              from "$targetDir/${task}_metrics/latent_mixing/main.nf"                params(params)
include  { ari }                        from "$targetDir/${task}_metrics/ari/main.nf"                          params(params)
include  { asw_batch }                  from "$targetDir/${task}_metrics/asw_batch/main.nf"                    params(params)
include  { asw_label }                  from "$targetDir/${task}_metrics/asw_label/main.nf"                    params(params)
include  { nmi }                        from "$targetDir/${task}_metrics/nmi/main.nf"                          params(params)
include  { cc_cons }                    from "$targetDir/${task}_metrics/cc_cons/main.nf"                      params(params)
include  { ti_cons }                    from "$targetDir/${task}_metrics/ti_cons/main.nf"                      params(params)
include  { graph_connectivity }         from "$targetDir/${task}_metrics/graph_connectivity/main.nf"           params(params)
include  { check_format }               from "$targetDir/${task}_metrics/check_format/main.nf"                 params(params)
include  { extract_scores }             from "$targetDir/common/extract_scores/main.nf"                        params(params)
include  { bind_tsv_rows }              from "$targetDir/common/bind_tsv_rows/main.nf"                         params(params)
include  { getDatasetId as get_id_predictions; getDatasetId as get_id_solutions } from "$srcDir/common/workflows/anndata_utils.nf"

params.datasets = "output/public_datasets/$task/**.h5ad"
workflow pilot_wf {
  main:

  // get input files for methods
  def inputs = 
    Channel.fromPath(params.datasets)
      | map { [ it.getParent().baseName, it ] }
      | filter { !it[1].name.contains("output_solution") }
      | groupTuple
      | map { id, datas -> 
        def fileMap = datas.collectEntries { [ (it.name.split(/\./)[-2].replace("output_", "input_")), it ]}
        [ id, fileMap, params ]
      }
  
  // get solutions
  def solution = 
    Channel.fromPath(params.datasets)
      | map { [ it.getParent().baseName, it ] }
      | filter { it[1].name.contains("output_solution") }

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
  def b3 = inputs 
    | baseline_totalvi
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_totalvi", [ input_prediction: pred, input_solution: sol ], params ]}
  def b4 = inputs 
    | baseline_mnn
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_mnn", [ input_prediction: pred, input_solution: sol ], params ]}
  def b5 = inputs 
    | baseline_newwave
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_newwave", [ input_prediction: pred, input_solution: sol ], params ]}

  def d0 = inputs 
    | dummy_random
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_dummy_random", [ input_prediction: pred, input_solution: sol ], params ]}
  def d1 = inputs 
    | dummy_zeros
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_dummy_zeros", [ input_prediction: pred, input_solution: sol ], params ]}
  def d2 = solution
    | map { id, input -> [ id, input, params ] }  
    | dummy_solution
    | join(solution)
    | map { id, pred, params, sol -> [ id + "_dummy_solution", [ input_prediction: pred, input_solution: sol ], params ]}

  def predictions = b0.mix(b1, b2, b3, b4, b5, d0, d1, d2)

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
    | (rfoob & latent_mixing & ari & asw_batch & asw_label & nmi & cc_cons & ti_cons & graph_connectivity & check_format)
    | mix
    | toList()
    | map{ [ it.collect{it[1]} ] }
    | combine(metricsMeta)
    | combine(solutionsMeta)
    | map{ [ "output", [ input: it[0], metric_meta: it[1], dataset_meta: it[2] ], params ] }
    | extract_scores
}
