nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"

// include  { generate_datasets }        from "$srcDir/common/workflows/generate_datasets/main.nf"                     params(params)
include  { prepare_task1_dataset }    from "$targetDir/predict_modality_datasets/prepare_task1_dataset/main.nf"     params(params)
include  { baseline_randomforest }    from "$targetDir/predict_modality_methods/baseline_randomforest/main.nf"      params(params)
include  { calculate_task1_metrics }  from "$targetDir/predict_modality_metrics/calculate_task1_metrics/main.nf"    params(params)


workflow run_task1_benchmark {
    main:
    output_ = Channel.fromPath("output/common_datasets/**.h5ad")
      | map { [ it.getParent().baseName, it, params ] }
      | view { "Received dataset '${it[0]}' at: ${it[1]}" }
      | prepare_task1_dataset
      | baseline_randomforest
      | calculate_task1_metrics
      
    emit: output_
}


// generate_datasets
