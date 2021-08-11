nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"

include  { baseline_dr_nn_knn }          from "$targetDir/match_modality_methods/baseline_dr_nn_knn/main.nf"          params(params)
include  { baseline_procrustes_knn }     from "$targetDir/match_modality_methods/baseline_procrustes_knn/main.nf"     params(params)
include  { dummy_constant }              from "$targetDir/match_modality_methods/dummy_constant/main.nf"              params(params)
include  { dummy_random }                from "$targetDir/match_modality_methods/dummy_random/main.nf"                params(params)
include  { calculate_auroc }             from "$targetDir/match_modality_metrics/calculate_auroc/main.nf"             params(params)
include  { extract_scores }              from "$targetDir/common/extract_scores/main.nf"                               params(params)

workflow pilot_wf {
  main:
  inputs = 
    Channel.fromPath("output/public_datasets/match_modality/**.h5ad")
      | map { [ it.getParent().baseName, it ] }
      | filter { !it[1].name.contains("output_solution") && !it[1].name.contains("output_test_sol") }
      | groupTuple
      | map { id, datas -> 
        def fileMap = datas.collectEntries { [ (it.name.split(/\./)[-2].replace("output_", "input_")), it ]}
        [ id, fileMap, params ]
      }
  solution = 
    Channel.fromPath("output/public_datasets/match_modality/**.h5ad")
      | map { [ it.getParent().baseName, it ] }
      | filter { it[1].name.contains("output_solution") || it[1].name.contains("output_test_sol") }

  // for now, code needs one of these code blocks per method.
  out0 = inputs 
    | baseline_dr_nn_knn
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_dr_nn_knn", [ input_prediction: pred, input_solution: sol ], params ]}

  out1 = inputs 
    | baseline_procrustes_knn
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_procrustes_knn", [ input_prediction: pred, input_solution: sol ], params ]}

  out2 = inputs 
    | dummy_constant
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_constant", [ input_prediction: pred, input_solution: sol ], params ]}

  out3 = inputs 
    | dummy_random
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_random", [ input_prediction: pred, input_solution: sol ], params ]}

  out0.mix(out1, out2, out3)
    | view{ [ "BASELINE", it[0], it[1] ] }
    | calculate_auroc
    | view{ [ "METRIC", it[0], it[1] ] }
    | map { it[1] }
    | toList()
    | map { [ "match_modality", it, params ] }
    | extract_scores
}
