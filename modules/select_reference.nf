process SELECT_REFERENCE {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/biocontainers/mulled-v2-77320db00eefbbf8c599692102c3d387a37ef02a:08144a66f00dc7684fad061f1466033c0176e7ad-0'
    
    input:
    tuple val(meta), path(db_covstats)

    output:
    tuple val(meta), path("*_refs.tsv"),            optional: true, emit: refs_tsv
    tuple val(meta), path("*_failed_assembly.tsv"), optional: true, emit: failed_assembly_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    select_reference.py \\
        -bbmap_covstats ${db_covstats} \\
        -b ${prefix} \\
        ${args}
    """
}
