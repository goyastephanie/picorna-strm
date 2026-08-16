//
// Soft-clip amplicon primer sequence from an alignment (ivar trim).
//
// Primer-derived bases reflect the primer, not the sample, so they are removed
// before consensus calling and before the per-position VCF. Trimming is done by
// COORDINATE (via the BED from MAKE_PRIMER_BED), not by matching the primer
// sequence in the read, so reads carrying only part of the primer are handled
// correctly too.
//
// iVar >= 1.4.1 is required: 1.4.1 fixed the -m parameter in trim, 1.4.2 a hang
// with -f, and 1.4.3 incorrect read positioning after trimming unpaired
// reverse reads.
//
process PRIMER_TRIM {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_medium'
    container 'quay.io/biocontainers/ivar:1.4.4--h077b44d_0'

    input:
    tuple val(meta), val(ref_info), path(bam), path(bai)
    tuple val(meta2), val(ref_info2), path(bed)

    output:
    tuple val(meta), val(ref_info), path("*.primertrim.bam"), path("*.primertrim.bam.bai"), emit: bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-e'
    def prefix = task.ext.prefix ?: "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    """
    if [ -s ${bed} ]; then
        # -e keeps reads that contain no primer, so coverage away from the
        # primer-binding site is preserved.
        ivar trim \\
            -i ${bam} \\
            -b ${bed} \\
            -p ${prefix}.trim \\
            ${args}

        samtools sort -@ ${task.cpus} -m 2G -o ${prefix}.primertrim.bam ${prefix}.trim.bam
    else
        echo "Empty primer BED (primer not found in this reference); passing the alignment through untrimmed."
        cp ${bam} ${prefix}.primertrim.bam
    fi

    samtools index -@ ${task.cpus} ${prefix}.primertrim.bam
    """
}
