nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"

include  { baseline_randomforest }       from "$targetDir/predict_modality_methods/baseline_randomforest/main.nf"      params(params)
include  { baseline_linearmodel }        from "$targetDir/predict_modality_methods/baseline_linearmodel/main.nf"       params(params)
include  { baseline_knearestneighbors }  from "$targetDir/predict_modality_methods/baseline_knearestneighbors/main.nf" params(params)
include  { dummy_zeros }                 from "$targetDir/predict_modality_methods/dummy_zeros/main.nf"                params(params)
include  { calculate_cor }               from "$targetDir/predict_modality_metrics/calculate_cor/main.nf"              params(params)
include  { extract_scores }              from "$targetDir/common/extract_scores/main.nf"                               params(params)

workflow pilot_wf {
  main:
  inputs = 
    Channel.fromPath("output/public_datasets/predict_modality/**.h5ad")
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
  rfOut = inputs 
    | map { id, data -> [ id, data, params ] }
    | baseline_randomforest
    | join(solution)
    | map { id, pred, params, sol -> [ id + "_rf", [ input_prediction: pred, input_solution: sol ], params ]}

  lmOut = inputs 
    | map { id, data -> [ id, data, params ] }
    | baseline_linearmodel
    | join(solution)
    | map { id, pred, params, sol -> [ id + "_lm", [ input_prediction: pred, input_solution: sol ], params ]}

  knnOut = inputs 
    | map { id, data -> [ id, data, params ] }
    | baseline_knearestneighbors
    | join(solution)
    | map { id, pred, params, sol -> [ id + "_knn", [ input_prediction: pred, input_solution: sol ], params ]}

  dzOut = inputs 
    | map { id, data -> [ id, data, params ] }
    | dummy_zeros
    | join(solution)
    | map { id, pred, params, sol -> [ id + "_dz", [ input_prediction: pred, input_solution: sol ], params ]}

  rfOut.mix(lmOut, knnOut, dzOut)
    // | view{ [ "BASELINE", it[0], it[1] ] }
    | calculate_cor
    // | view{ [ "METRIC", it[0], it[1] ] }
    | map { it[1] }
    | toList()
    | map { [ "predict_modality", it, params ] }
    | extract_scores
}
