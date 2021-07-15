nextflow.enable.dsl=2

targetDir = "${params.rootDir}/target/nextflow"

include  { prepare_task1_dataset }       from "$targetDir/predict_modality_datasets/prepare_task1_dataset/main.nf"     params(params)

// params.prepare_task1_dataset__max_mod1_columns = 1000
params.prepare_task1_dataset__max_mod2_columns = 1000

workflow {
  main:
  modsAndSolutions = 
    Channel.fromPath("output/common_datasets/**.h5ad")
      | map { [ it.getParent().baseName, it ] }
      | groupTuple
      | flatMap { id, datas -> 
        def fileMap = datas.collectEntries { [ (it.name.split(/\./)[-2]), it ]}
        [
          [ id + "_rna", [ input_mod1: fileMap.output_rna, input_mod2: fileMap.output_mod2 ], params ],
          [ id + "_mod2", [ input_mod1: fileMap.output_mod2, input_mod2: fileMap.output_rna ], params ] 
        ]
      }
      | prepare_task1_dataset
}
