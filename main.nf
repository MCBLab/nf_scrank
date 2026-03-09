#!/usr/bin/env nextflow

params.seurat_obj = "vehicle_final_server.rds"
params.targets_file = "alvos.txt"
params.species = "mouse"
params.outdir = "./resultados"

seurat_ch = Channel.fromPath(params.seurat_obj)
targets_ch = Channel.fromPath(params.targets_file)

process split_seurat {
    publishDir "${params.outdir}/inputs", mode: 'copy'
    input: path seurat_file; path targets_file
    output: path "nextflow_inputs/*.rds"
    script: "Rscript ${projectDir}/01_split_seurat.R $seurat_file $targets_file ${params.species}"
}

process run_genie3 {
    cpus 4
    input: path split_rds
    output: path "*_weights.rds"
    script: "Rscript ${projectDir}/02_run_genie3_worker.R $split_rds 4"
}

process merge_scores {
    publishDir "${params.outdir}/final", mode: 'copy'
    input: path weight_files
    output: path "*.csv"
    script: "Rscript ${projectDir}/03_merge_scores.R 'Smo' $weight_files"
}

workflow {
    split_parts = split_seurat(seurat_ch, targets_ch)
    weights = run_genie3(split_parts.flatten())
    merge_scores(weights.collect())
}
