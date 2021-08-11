nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"

include  { baseline_lmds }              from "$targetDir/joint_embedding_methods/baseline_lmds/main.nf"                params(params)
include  { baseline_pca }               from "$targetDir/joint_embedding_methods/baseline_pca/main.nf"                 params(params)
include  { baseline_umap }              from "$targetDir/joint_embedding_methods/baseline_umap/main.nf"                params(params)
include  { dummy_random }               from "$targetDir/joint_embedding_methods/dummy_random/main.nf"                 params(params)
include  { dummy_single_mod }           from "$targetDir/joint_embedding_methods/dummy_single_mod/main.nf"             params(params)
include  { dummy_zeros }                from "$targetDir/joint_embedding_methods/dummy_zeros/main.nf"                  params(params)
include  { totalvi }                    from "$targetDir/joint_embedding_methods/totalvi/main.nf"                      params(params)
include  { calculate_rf_oob }           from "$targetDir/joint_embedding_metrics/calculate_rf_oob/main.nf"             params(params)
include  { calculate_totalVI_metrics }  from "$targetDir/joint_embedding_metrics/calculate_totalVI_metrics/main.nf"    params(params)
include  { ari }                        from "$targetDir/joint_embedding_metrics/ari/main.nf"                          params(params)
include  { asw_batch }                  from "$targetDir/joint_embedding_metrics/asw_batch/main.nf"                    params(params)
include  { asw_label }                  from "$targetDir/joint_embedding_metrics/asw_label/main.nf"                    params(params)
include  { nmi }                        from "$targetDir/joint_embedding_metrics/nmi/main.nf"                          params(params)
include  { extract_scores }             from "$targetDir/common/extract_scores/main.nf"                                params(params)

workflow pilot_wf {
  main:
  inputs = 
    Channel.fromPath("output/public_datasets/joint_embedding/**.h5ad")
      | map { [ it.getParent().baseName, it ] }
      | filter { !it[1].name.contains("output_solution") && !it[1].name.contains("output_test_sol") }
      | groupTuple
      | map { id, datas -> 
        def fileMap = datas.collectEntries { [ (it.name.split(/\./)[-2].replace("output_", "input_")), it ]}
        [ id, fileMap, params ]
      }
  solution = 
    Channel.fromPath("output/public_datasets/joint_embedding/**.h5ad")
      | map { [ it.getParent().baseName, it ] }
      | filter { it[1].name.contains("output_solution") || it[1].name.contains("output_test_sol") }

  // for now, code needs one of these code blocks per method.
  b0 = inputs 
    | baseline_lmds
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_lmds", [ input_prediction: pred, input_solution: sol ], params ]}
  b1 = inputs 
    | baseline_pca
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_pca", [ input_prediction: pred, input_solution: sol ], params ]}
  b2 = inputs 
    | baseline_umap
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_baseline_umap", [ input_prediction: pred, input_solution: sol ], params ]}

  d0 = inputs 
    | dummy_random
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_dummy_random", [ input_prediction: pred, input_solution: sol ], params ]}
  d1 = inputs 
    | dummy_single_mod
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_dummy_single_mod", [ input_prediction: pred, input_solution: sol ], params ]}
  d2 = inputs 
    | dummy_zeros
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_dummy_zeros", [ input_prediction: pred, input_solution: sol ], params ]}

  m0 = inputs 
    | totalvi
    | join(solution) 
    | map { id, pred, params, sol -> [ id + "_totalvi", [ input_prediction: pred, input_solution: sol ], params ]}

  b0.mix(b1, b2, d0, d1, d2, m0)
    | view{ [ "BASELINE", it[0], it[1] ] }
    | (calculate_rf_oob & calculate_totalVI_metrics & ari & asw_batch & asw_label & nmi)
    | mix
    | view{ [ "METRIC", it[0], it[1] ] }
    | map { it[1] }
    | toList()
    | map { [ "joint_embedding", it, params ] }
    | extract_scores
}
