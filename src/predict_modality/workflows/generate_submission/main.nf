nextflow.enable.dsl=2

include { method } from "$launchDir/target/nextflow/main.nf" params(params)

params.datasets = "s3://neurips2021-multimodal-public-datasets/predict_modality/**.h5ad"

workflow {
  main:
  print(params.datasets)
  Channel.fromPath(params.datasets)
    | map { [ it.getParent().baseName, it ] }
    | filter { !it[1].name.contains("output_test_mod2") }
    | groupTuple
    | map { id, datas -> 
      def fileMap = datas.collectEntries { [ (it.name.split(/\./)[-2].replace("output_", "input_")), it ]}
      [ id, fileMap, params ]
    }
    | method
}
