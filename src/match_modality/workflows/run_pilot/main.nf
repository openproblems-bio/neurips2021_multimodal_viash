nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"
task = "match_modality"

// always fails
// include  { baseline_babel_knn }          from "$targetDir/${task}_methods/baseline_babel_knn/main.nf"          params(params)
include  { baseline_dr_knnr_knn }        from "$targetDir/${task}_methods/baseline_dr_knnr_knn/main.nf"        params(params)
include  { baseline_dr_knnr_cbf }        from "$targetDir/${task}_methods/baseline_dr_knnr_cbf/main.nf"        params(params)
include  { baseline_newwave_knnr_cbf }   from "$targetDir/${task}_methods/baseline_newwave_knnr_cbf/main.nf"   params(params)
include  { baseline_newwave_knnr_knn }   from "$targetDir/${task}_methods/baseline_newwave_knnr_knn/main.nf"   params(params)
include  { baseline_procrustes_knn }     from "$targetDir/${task}_methods/baseline_procrustes_knn/main.nf"     params(params)
include  { baseline_linear_knn }         from "$targetDir/${task}_methods/baseline_linear_knn/main.nf"         params(params)
include  { dummy_constant }              from "$targetDir/${task}_methods/dummy_constant/main.nf"              params(params)
include  { dummy_random }                from "$targetDir/${task}_methods/dummy_random/main.nf"                params(params)
include  { dummy_solution }              from "$targetDir/${task}_methods/dummy_solution/main.nf"              params(params)
include  { dummy_semisolution }          from "$targetDir/${task}_methods/dummy_semisolution/main.nf"          params(params)
include  { aupr }                        from "$targetDir/${task}_metrics/aupr/main.nf"                        params(params)
include  { match_probability }           from "$targetDir/${task}_metrics/match_probability/main.nf"           params(params)
include  { check_format }                from "$targetDir/${task}_metrics/check_format/main.nf"                params(params)
include  { final_scores }                from "$targetDir/${task}_results/final_scores/main.nf"                params(params)
include  { bind_tsv_rows }               from "$targetDir/common/bind_tsv_rows/main.nf"                        params(params)

params.datasets = "output/public_datasets/$task/**.h5ad"
params.meta_datasets = "${params.rootDir}/results/meta_datasets.tsv"

workflow pilot_wf {
  main:

  // get input files for methods
  def inputs = 
    Channel.fromPath(params.datasets)
      | map { [ it.getParent().baseName, it ] }
      | filter { !it[1].name.contains("output_test_sol") }
      | groupTuple
      | map { id, datas -> 
        def fileMap = datas.collectEntries { [ (it.name.split(/\./)[-2].replace("output_", "input_")), it ]}
        [ id, fileMap, params ]
      }
  
  // get solutions
  def solution = 
    Channel.fromPath(params.datasets)
      | map { [ it.getParent().baseName, it ] }
      | filter { it[1].name.contains("output_test_sol") }

  // for now, code needs one of these code blocks per method.
  def b0 = inputs 
    | baseline_procrustes_knn
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_procrustes_knn", [ input_prediction: pred, input_solution: sol ], params ]}
  def b1 = inputs 
    | baseline_dr_knnr_knn
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_dr_knnr_knn", [ input_prediction: pred, input_solution: sol ], params ]}
  // def b2 = inputs 
  //   | baseline_newwave_knnr_cbf
  //   | join(solution) 
  //   | map { id, pred, params, sol -> [ id + "_baseline_newwave_knnr_cbf", [ input_prediction: pred, input_solution: sol ], params ]}
  // def b3 = inputs 
  //   | baseline_newwave_knnr_knn
  //   | join(solution) 
  //   | map { id, pred, params, sol -> [ id + "_baseline_newwave_knnr_knn", [ input_prediction: pred, input_solution: sol ], params ]}
  // def b4 = inputs 
  //   | baseline_babel_knn
  //   | join(solution) 
  //   | map { id, pred, params, sol -> [ id + "_baseline_babel_knn", [ input_prediction: pred, input_solution: sol ], params ]}
  def b5 = inputs 
    | baseline_linear_knn
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_linear_knn", [ input_prediction: pred, input_solution: sol ], params ]}
  def b6 = inputs 
    | baseline_dr_knnr_cbf
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_dr_knnr_cbf", [ input_prediction: pred, input_solution: sol ], params ]}

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
    | dummy_solution
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_dummy_solution", [ input_prediction: pred, input_solution: sol ], params ]}
  def d3 = solution
    | map { id, input -> [ id, input, params ] } 
    | dummy_semisolution
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_dummy_semisolution", [ input_prediction: pred, input_solution: sol ], params ]}

  def predictions = b0.mix(b1, b5, b6, d0, d1, d2, d3)

  // fetch dataset ids in predictions and in solutions
  def datasetsMeta = 
    Channel.fromPath(params.meta_datasets)
  
  // create metrics meta
  def metricsMeta = 
    Channel.fromPath("$srcDir/$task/**/metric_meta*.tsv")
      | toList()
      | map{ [ "meta", it, params ] }
      | bind_tsv_rows
      | map{ it[1] }

  // compute metrics & combine results
  predictions
    | (aupr & match_probability & check_format)
    | mix
    | toList()
    | map{ [ it.collect{it[1]} ] }
    | combine(metricsMeta)
    | combine(datasetsMeta)
    | map{ [ "output", [ input: it[0], metric_meta: it[1], dataset_meta: it[2] ], params ] }
    | final_scores
}
