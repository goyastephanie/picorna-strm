//
// Locate the amplicon primer in a reference and emit a BED for ivar trim.
//
// Split out from PRIMER_TRIM so that the trimming step can run in the clean
// iVar biocontainer (iVar + samtools, no python) while the search runs in a
// python image. The primer sits at a different coordinate in every genotype,
// and at a different coordinate again in the sample's own first consensus, so
// this is recomputed per alignment rather than supplied by the user.
//
process MAKE_PRIMER_BED {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_single'
    container 'quay.io/biocontainers/python:3.8.3'

    input:
    tuple val(meta), val(ref_info), path(ref)
    val primer_fwd

    output:
    tuple val(meta), val(ref_info), path("*_primers.bed"), emit: bed

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    """
    make_primer_bed.py \\
        --ref ${ref} \\
        --primer ${primer_fwd} \\
        --max_mismatch ${params.primer_max_mismatch} \\
        --out ${prefix}_primers.bed
    """
}
