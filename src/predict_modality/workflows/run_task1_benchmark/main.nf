nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"

include  { baseline_randomforest }       from "$targetDir/predict_modality_methods/baseline_randomforest/main.nf"      params(params)
include  { baseline_linearmodel }        from "$targetDir/predict_modality_methods/baseline_linearmodel/main.nf"       params(params)
include  { baseline_knearestneighbors }  from "$targetDir/predict_modality_methods/baseline_knearestneighbors/main.nf" params(params)
include  { calculate_task1_metrics }     from "$targetDir/predict_modality_metrics/calculate_task1_metrics/main.nf"    params(params)
include  { extract_scores }              from "$targetDir/common/extract_scores/main.nf"                               params(params)

workflow run_task1_benchmark {
  main:
  modsAndSolutions = 
    Channel.fromPath("output/task1_datasets/**.h5ad")
      | map { [ it.getParent().baseName, it ] }
      | groupTuple
      | map { id, datas -> 
        def fileMap = datas.collectEntries { [ (it.name.split(/\./)[-2]), it ]}
        [ id, fileMap ]
      }

  // for now, code needs one of these code blocks per method.
  rfOut = modsAndSolutions 
    | map { id, data -> [ id, [ input_mod1: data.output_mod1, input_mod2: data.output_mod2 ], params ] }
    | baseline_randomforest
    | join(modsAndSolutions) 
    | map { id, pred, params, modsAndSols -> [ id + "_rf", [ input_prediction: pred, input_solution: modsAndSols.output_solution ], params ]}

  lmOut = modsAndSolutions 
    | map { id, data -> [ id, [ input_mod1: data.output_mod1, input_mod2: data.output_mod2 ], params ] }
    | baseline_linearmodel
    | join(modsAndSolutions) 
    | map { id, pred, params, modsAndSols -> [ id + "_lm", [ input_prediction: pred, input_solution: modsAndSols.output_solution ], params ]}

  knnOut = modsAndSolutions 
    | map { id, data -> [ id, [ input_mod1: data.output_mod1, input_mod2: data.output_mod2 ], params ] }
    | baseline_knearestneighbors
    | join(modsAndSolutions) 
    | map { id, pred, params, modsAndSols -> [ id + "_knn", [ input_prediction: pred, input_solution: modsAndSols.output_solution ], params ]}

  rfOut.mix(lmOut, knnOut)
    // | view{ [ "BASELINE", it[0], it[1] ] }
    | calculate_task1_metrics
    // | view{ [ "METRIC", it[0], it[1] ] }
    | map { it[1] }
    | toList()
    | map { [ "task1", it, params ] }
    | extract_scores
}
