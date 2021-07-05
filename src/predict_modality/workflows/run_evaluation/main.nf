nextflow.enable.dsl=2

workflowDir = "${params.rootDir}/src/utils"
targetDir = "${params.rootDir}/target/nextflow"

include  { download_10x_dataset }     from "$targetDir/common_datasets/download_10x_dataset/main.nf"                params(params)
include  { dyngen }                   from "$targetDir/common_datasets/dyngen/main.nf"                              params(params)
// include  { prepare_task1_dataset }    from targetDir + "/predict_modality_datasets/prepare_task1_dataset/main.nf"     params(params)
// include  { baseline_randomforest }    from targetDir + "/predict_modality_methods/baseline_randomforest/main.nf"      params(params)
// include  { calculate_task1_metrics }  from targetDir + "/predict_modality_metrics/calculate_task1_metrics/main.nf"    params(params)
// include  { overrideOptionValue }      from  workflowDir + "/utils/utils.nf"                        params(params)
include  { overrideOptionValue }     from  "${params.rootDir}/src/utils/workflows/utils.nf"


/* BD Rhapsody WTA - multi-sample workflow
 * 
 * consumed params:
 *   tsv:                           a tsv for processing multiple input files. tsv must contain two columns, 'id' and 'input'. 
 *                                  'id' is a sample id for one or more fastq files. 
 *                                  'input' is one or more fastq paths, separated with semicolons, paths may be globs
 *   reference_genome:              a path to STAR index as a tar.gz file
 *   transcriptome_annotation:      a path to GTF annotation file
 *   output                         a publish dir for the output h5ad files
 * output format:               [ id, h5ad, params ]
 *   value id:                      a sample id for one or more fastq files
 *   value h5ad:                    h5ad object of mapped fastq reads
 *   value params:                  the params object, which may already have sample specific overrides
 * publishes:
 *   the output h5ad files
 */
workflow dyngen_datasets {
    main:
    if (!params.containsKey("tsv") || params.tsv == "") {
        exit 1, "ERROR: Please provide a --tsv parameter. The tsv must include two columns 'id' and 'input'."
    }
    
    output_ = Channel.fromPath(file(params.tsv)) \
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

workflow public_10x_datasets {
    main:
    if (!params.containsKey("tsv") || params.tsv == "") {
        exit 1, "ERROR: Please provide a --tsv parameter. The tsv must include two columns 'id' and 'input'."
    }
    
    output_ = Channel.fromPath(file(params.tsv)) \
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
