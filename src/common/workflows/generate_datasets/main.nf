nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"

include  { download_10x_dataset }   from "$targetDir/common_datasets/download_10x_dataset/main.nf"                params(params)
include  { dyngen }                 from "$targetDir/common_datasets/simulate_dyngen_dataset/main.nf"                              params(params)
include  { overrideOptionValue }    from  "$srcDir/common/workflows/utils.nf"


workflow generate_dyngen_datasets {
    main:
    output_ = Channel.fromPath(file("$srcDir/common/datasets/simulate_dyngen_dataset/input.tsv")) \
        | splitCsv(header: true, sep: "\t") \
        | map { tsv -> [ tsv.id, tsv, params ] } \
        | map { overrideOptionValue(it, "dyngen", "id", it[1].id) } \
        | map { overrideOptionValue(it, "dyngen", "backbone", it[1].backbone) } \
        | map { overrideOptionValue(it, "dyngen", "num_cells", it[1].num_cells) } \
        | map { overrideOptionValue(it, "dyngen", "num_genes", it[1].num_genes) } \
        | map { overrideOptionValue(it, "dyngen", "num_simulations", it[1].num_simulations) } \
        | map { overrideOptionValue(it, "dyngen", "store_chromatin", it[1].store_chromatin) } \
        | map { overrideOptionValue(it, "dyngen", "store_protein", it[1].store_protein) } \
        | map { [ it[0], [], it[2] ] } \
        | filter{ it[0] ==~ /.*_small/ } \
        | view{ [ "DEBUG2", it[0], it[1] ] } \
        | dyngen \
        | view{ [ "DEBUG3", it[0], it[1] ] }

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
        | view{ [ "DEBUG2", it[0], it[1] ] } \
        | download_10x_dataset \
        | view{ [ "DEBUG3", it[0], it[1] ] }

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
      | mix
      
    emit: output_
}
