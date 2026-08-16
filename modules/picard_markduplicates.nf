process PICARD_MARKDUPLICATES {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_medium'
    container 'quay.io/biocontainers/picard:3.0.0--hdfd78af_1'

    input:
    tuple val(meta), val(ref_info), path(bam), path(bai)

    output:
    tuple val(meta), val(ref_info), path("*.markdup.bam"), path("*.markdup.bam.bai"), emit: bam
    tuple val(meta), val(ref_info), path("*.metrics.txt"), optional: true,            emit: metrics

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    picard \\
        -Xmx${avail_mem}M \\
        MarkDuplicates \\
        --REMOVE_DUPLICATES true \\
        --CREATE_INDEX true \\
        --ASSUME_SORT_ORDER coordinate \\
        $args \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.markdup.bam \\
        --METRICS_FILE ${prefix}.MarkDuplicates.metrics.txt

    # Picard writes the index as "<prefix>.markdup.bai"; rename it to the
    # "<bam>.bai" convention used everywhere else in the pipeline so samtools
    # picks it up unambiguously downstream.
    if [ -f ${prefix}.markdup.bai ]; then
        mv ${prefix}.markdup.bai ${prefix}.markdup.bam.bai
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    """
    touch ${prefix}.markdup.bam
    touch ${prefix}.markdup.bam.bai
    touch ${prefix}.MarkDuplicates.metrics.txt
    """
}
