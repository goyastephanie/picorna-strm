//
// Per-position pileup VCF against the sample's own consensus.
//
// This is NOT a variant caller: bcftools mpileup is run WITHOUT bcftools call,
// so every position is reported with its allelic depths (FORMAT/AD) and total
// depth (FORMAT/DP). That makes the file a complete record of intra-host
// (intra-genotype) diversity, from which allele frequencies can be computed at
// any site downstream -- including sites no caller would have flagged.
//
// The alignment is the one the FINAL consensus was built from (reads vs the
// first consensus), already primer-trimmed and/or deduplicated according to
// --mode, and the reference is that same first consensus.
//
process BCFTOOLS_MPILEUP {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_medium'
    container 'quay.io/biocontainers/bcftools:1.17--haef29d1_0'

    input:
    tuple val(meta), val(ref_info), path(bam), path(bai)
    tuple val(meta2), val(ref_info2), path(ref)

    output:
    tuple val(meta), val(ref_info), path("*.vcf"), optional: true, emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: "-a FORMAT/AD,FORMAT/DP --max-depth ${params.vcf_max_depth} --min-BQ ${params.vcf_min_bq} --min-MQ ${params.vcf_min_mq}"
    def prefix = task.ext.prefix ?: "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    """
    if [[ \$(basename "$bam") = "FAILED.sorted.bam" ]]; then
        echo "Skipping pileup for ${prefix}; alignment failed thresholds previously"
        exit 0
    fi

    # Work on a local copy: bcftools mpileup builds the .fai next to the
    # reference, and the staged file is read-only. (The bcftools container has
    # no samtools, and bcftools has no 'faidx' subcommand.)
    cp ${ref} ref_for_pileup.fa

    bcftools mpileup \\
        -f ref_for_pileup.fa \\
        ${args} \\
        -Ov \\
        -o ${prefix}.vcf \\
        ${bam}
    """
}
