nextflow.enable.dsl=2

srcDir = "${params.rootDir}/src"
targetDir = "${params.rootDir}/target/nextflow"

include  { download_10x_dataset }               from "$targetDir/common_datasets/download_10x_dataset/main.nf"                params(params)
include  { simulate_dyngen_dataset }            from "$targetDir/common_datasets/simulate_dyngen_dataset/main.nf"             params(params)
include  { download_azimuth_dataset }           from "$targetDir/common_datasets/download_azimuth_dataset/main.nf"            params(params)
include  { download_totalvi_spleen_lymph }      from "$targetDir/common_datasets/download_totalvi_spleen_lymph/main.nf"       params(params)
include  { download_totalvi_10x_datasets }      from "$targetDir/common_datasets/download_totalvi_10x_datasets/main.nf"       params(params)
include  { quality_control }                    from "$targetDir/common_process_dataset/quality_control/main.nf"              params(params)
include  { split_traintest }                    from "$targetDir/common_process_dataset/split_traintest/main.nf"              params(params)
include  { simulate_batch }                     from "$targetDir/common_process_dataset/simulate_batch/main.nf"               params(params)
include  { pseudotime_order }                   from "$targetDir/common_process_dataset/pseudotime_order/main.nf"             params(params)
include  { overrideOptionValue }                from "$srcDir/common/workflows/utils.nf"

def flattenMap(entry) {
  res = [:]
  entry.each{it.each{ res[it.key] = it.value }}
  return res
}

workflow generate_dyngen_datasets {
    main:
    def cacheDir = file("${workflow.homeDir}/.local/share/R/dyngen")
    output_ = Channel.fromPath(file("$srcDir/common/datasets/simulate_dyngen_dataset/input.tsv"))
        | splitCsv(header: true, sep: "\t")
        | map { tsv -> [ tsv.id, tsv, params ] }
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "id", it[1].id) }
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "num_cells", it[1].num_cells) }
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "num_genes", it[1].num_genes) }
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "num_simulations", it[1].num_simulations) }
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "store_chromatin", it[1].store_chromatin) }
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "store_protein", it[1].store_protein) }
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "num_threads", it[1].num_threads) }
        | map { overrideOptionValue(it, "simulate_dyngen_dataset", "seed", it[1].seed) }
        | map { [ it[0], [cache_dir: cacheDir], it[2] ] }
        | simulate_dyngen_dataset

    emit: output_
}

workflow generate_public_10x_datasets {
    main:
    output_ = Channel.fromPath(file("$srcDir/common/datasets/download_10x_dataset/input.tsv"))
        | splitCsv(header: true, sep: "\t")
        | map { tsv -> [ tsv.id, tsv, params ] }
        | map { overrideOptionValue(it, "download_10x_dataset", "id", it[1].id) }
        | map { overrideOptionValue(it, "download_10x_dataset", "input", it[1].input) }
        | map { overrideOptionValue(it, "download_10x_dataset", "organism", it[1].organism) }
        | map { [ it[0], [], it[2] ] }
        | download_10x_dataset

    emit: output_
}

workflow generate_azimuth_datasets {
    main:
    output_ = Channel.fromPath(file("$srcDir/common/datasets/download_azimuth_dataset/input.tsv"))
        | splitCsv(header: true, sep: "\t")
        | map { tsv -> [ tsv.id, tsv, params ] }
        | map { overrideOptionValue(it, "download_azimuth_dataset", "id", it[1].id) }
        | map { overrideOptionValue(it, "download_azimuth_dataset", "input_count", it[1].input_count) }
        | map { overrideOptionValue(it, "download_azimuth_dataset", "input_meta", it[1].input_meta) }
        | map { overrideOptionValue(it, "download_azimuth_dataset", "organism", it[1].organism) }
        | map { [ it[0], [], it[2] ] }
        | download_azimuth_dataset

    emit: output_
}

workflow generate_totalvi_spleen_lymph {
    main:
    output_ = Channel.fromPath(file("$srcDir/common/datasets/download_totalvi_spleen_lymph/input.tsv"))
        | splitCsv(header: true, sep: "\t")
        | map { tsv -> [ tsv.id, tsv, params ] }
        | map { overrideOptionValue(it, "download_totalvi_spleen_lymph", "id", it[1].id) }
        | map { overrideOptionValue(it, "download_totalvi_spleen_lymph", "input", it[1].input) }
        | map { overrideOptionValue(it, "download_totalvi_spleen_lymph", "organism", it[1].organism) }
        | map { [ it[0], [], it[2] ] }
        | download_totalvi_spleen_lymph

    emit: output_
}

workflow generate_totalvi_10x_datasets {
    main:
    output_ = Channel.fromPath(file("$srcDir/common/datasets/download_totalvi_10x_datasets/input.tsv"))
        | splitCsv(header: true, sep: "\t")
        | map { tsv -> [ tsv.id, tsv, params ] }
        | map { overrideOptionValue(it, "download_totalvi_10x_datasets", "id", it[1].id) }
        | map { overrideOptionValue(it, "download_totalvi_10x_datasets", "input", it[1].input) }
        | map { overrideOptionValue(it, "download_totalvi_10x_datasets", "organism", it[1].organism) }
        | map { [ it[0], [], it[2] ] }
        | download_totalvi_10x_datasets

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
    output_ = (generate_dyngen_datasets & generate_public_10x_datasets & generate_azimuth_datasets & generate_totalvi_spleen_lymph & generate_totalvi_10x_datasets)
      | mix
      | map { id, data, prms -> [ id, [ input_rna: data.output_rna, input_mod2: data.output_mod2 ], prms ]}
      // | view { ["DEBUG0", it[0], it[1] ] }
      | quality_control
      | map { id, data, prms -> [ id, [ input_rna: data.output_rna, input_mod2: data.output_mod2 ], prms ]}
      // | view { ["DEBUG1", it[0], it[1] ] }
      | pseudotime_order
      | map { id, data, prms -> [ id, [ input_rna: data.output_rna, input_mod2: data.output_mod2 ], prms ]}
      // | view { ["DEBUG1", it[0], it[1] ] }
      | simulate_batch
      | map { id, data, prms -> [ id, [ input_rna: data.output_rna, input_mod2: data.output_mod2 ], prms ]}
      // | view { ["DEBUG1", it[0], it[1] ] }
      | split_traintest
      | view { "Publishing dataset with ${it[0]} from ${it[1]}" }
      
    emit: output_
}
