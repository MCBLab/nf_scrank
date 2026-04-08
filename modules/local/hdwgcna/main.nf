process HDWGCNA {
    tag "${rds_file.baseName}"
    publishDir "${params.outdir}/hdwgcna_networks", mode: 'copy'
    
    container 'docker://leoshow21/hdwgcna:v2'
    
    queue 'amd-512'
    cpus params.n_cores
    memory '150 GB'
    
    input:
    path rds_file
    val n_cores
    
    output:
    path "*_weight_hdWGCNA_*.rds", emit: rank_obj, optional: true
    
    script:
    """
    hdwgcna.R ${rds_file} ${n_cores} ${params.target}
    """
}
