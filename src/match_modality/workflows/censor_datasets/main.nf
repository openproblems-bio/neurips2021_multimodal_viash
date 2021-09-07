nextflow.enable.dsl=2

targetDir = "${params.rootDir}/target/nextflow"

include  { censor_dataset } from "$targetDir/match_modality_datasets/censor_dataset/main.nf" params(params)


workflow {
  main:
  Channel.fromPath(params.datasets)
    | map { [ it.getParent().baseName, it ] }
    | groupTuple
    | flatMap { id, datas -> 
      def fileMap = datas.collectEntries { [ (it.name.split(/\./)[-2]), it ]}
      [
        [ id + "_rna", [ input_mod1: fileMap.output_rna, input_mod2: fileMap.output_mod2 ], params ],
        [ id + "_mod2", [ input_mod1: fileMap.output_mod2, input_mod2: fileMap.output_rna ], params ] 
      ]
    }
    | censor_dataset
}
