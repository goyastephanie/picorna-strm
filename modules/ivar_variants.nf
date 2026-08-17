//
// Variant calling against the sample's own final consensus (iVar variants).
//
// Runs on the same alignment used to build the final consensus, so the
// variants are primer-trimmed exactly like the consensus is.
// Emits iVar's native TSV with per-allele depth and frequency. The per-position
// VCF is produced separately by BCFTOOLS_MPILEUP.
//
process IVAR_VARIANTS {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    // samtools mpileup | ivar variants: both single-threaded
    label 'process_single'
    container 'quay.io/biocontainers/ivar:1.4.4--h077b44d_0'

    input:
    tuple val(meta), val(ref_info), path(bam), path(bai)
    tuple val(meta2), val(ref_info2), path(ref)

    output:
    tuple val(meta), val(ref_info), path("*.variants.tsv"), optional: true, emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "-q ${params.ivar_var_q} -t ${params.ivar_var_t} -m ${params.ivar_var_m}"
    def args2 = task.ext.args2 ?: "--count-orphans --no-BAQ --max-depth ${params.mpileup_max_depth} --min-BQ 0 -aa"
    def prefix = task.ext.prefix ?: "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    """
    if [[ \$(basename "$bam") = "FAILED.sorted.bam" ]]; then
        echo "Skipping variant calling for ${prefix}; alignment failed thresholds previously"
        exit 0
    fi

    samtools \\
        mpileup \\
        --reference ${ref} \\
        ${args2} \\
        ${bam} \\
        | ivar \\
            variants \\
            ${args} \\
            -r ${ref} \\
            -p ${prefix}.variants
    """
}
