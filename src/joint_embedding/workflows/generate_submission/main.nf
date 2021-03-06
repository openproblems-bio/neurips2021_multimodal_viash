nextflow.enable.dsl=2

include { method } from "$launchDir/target/nextflow/main.nf" params(params)

params.datasets = "s3://neurips2021-multimodal-public-datasets/joint_embedding/**.h5ad"

workflow {
  main:
  Channel.fromPath(params.datasets)
    | map { [ it.getParent().baseName, it ] }
    | filter { !it[1].name.contains("output_solution") }
    // | view { [ "DEBUG0", it[0], it[1] ]}
    | groupTuple
    | map { id, datas -> 
      def fileMap = datas.collectEntries { [ (it.name.split(/\./)[-2].replace("output_", "input_")), it ]}
      [ id, fileMap, params ]
    }
    | method
}