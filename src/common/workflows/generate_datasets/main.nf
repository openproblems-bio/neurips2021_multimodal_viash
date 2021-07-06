nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"

include  { download_10x_dataset }    from "$targetDir/common_datasets/download_10x_dataset/main.nf"                params(params)
include  { simulate_dyngen_dataset } from "$targetDir/common_datasets/simulate_dyngen_dataset/main.nf"             params(params)
include  { normalize }               from "$targetDir/common_datasets/normalize/main.nf"                           params(params)
include  { overrideOptionValue }     from "$srcDir/common/workflows/utils.nf"

workflow generate_dyngen_datasets {
    main:
    output_ = Channel.fromPath(file("$srcDir/common/datasets/simulate_dyngen_dataset/input.tsv")) \
        | splitCsv(header: true, sep: "\t") \
        | map { tsv -> [ tsv.id, tsv, params ] } \
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "id", it[1].id) } \
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "backbone", it[1].backbone) } \
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "num_cells", it[1].num_cells) } \
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "num_genes", it[1].num_genes) } \
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "num_simulations", it[1].num_simulations) } \
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "store_chromatin", it[1].store_chromatin) } \
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "store_protein", it[1].store_protein) } \
        | map { [ it[0], [], it[2] ] } \
        | filter { it[0] ==~ /.*_small/ } \
        | simulate_dyngen_dataset \
        | filter { it[1].containsKey("output") } \
        | map { [ it[0], it[1].output, it[2] ] }

    emit: output_
}

workflow generate_public_10x_datasets {
    main:
    output_ = Channel.fromPath(file("$srcDir/common/datasets/download_10x_dataset/input.tsv")) \
        | splitCsv(header: true, sep: "\t") \
        | map { tsv -> [ tsv.id, tsv, params ] } \
        | map { overrideOptionValue(it, "download_10x_dataset", "id", it[1].id) } \
        | map { overrideOptionValue(it, "download_10x_dataset", "input", it[1].input) } \
        | map { [ it[0], [], it[2] ] } \
        | download_10x_dataset

    emit: output_
}


/* data gen workflow
 * 
 * consumed params:
 *   <none>
 * output format:               [ id, h5ad, params ]
 *   value id:                      a dataset id
 *   value h5ad:                    the dataset h5ad file
 *   value params:                  the (dataset-specific) parameters
 * publishes:
 *   the output h5ad files
 */
workflow generate_datasets {
    main:
    output_ = (generate_dyngen_datasets & generate_public_10x_datasets) \
      | mix \
      | normalize \
      | view { "Publishing dataset with ${it[0]} from ${it[1]}" }
      
    emit: output_
}
