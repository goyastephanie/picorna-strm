//
// Locate the amplicon primer in a reference and emit a BED for ivar trim.
//
// Split out from PRIMER_TRIM so that the trimming step can run in the clean
// iVar biocontainer (iVar + samtools, no python) while the search runs in a
// python image.
//
// The search is repeated for each alignment pass because the two references
// are not the same sequence:
//
//   pass 1 (vs the database reference) -- the primer sits at its genomic
//     coordinate, which differs between genotypes (~pos 300-450 in the 5'UTR
//     for rhinovirus). It is found and trimmed.
//
//   pass 2 (vs the sample's first consensus) -- acts as a safety net. When
//     pass 1 worked, the primer region has zero depth, iVar writes Ns there and
//     IVAR_CONSENSUS strips leading Ns, so the first consensus begins just
//     downstream of the primer: the search returns an empty BED and nothing is
//     trimmed (correct -- reads still carrying the primer have those bases
//     hanging off the 5' end and the aligner soft-clips them anyway). If any
//     primer DID survive into the first consensus, it is caught here. Because
//     quality trimming leaves a variable amount of primer on each read, that
//     residual primer is often TRUNCATED at the very start of the consensus, so
//     the search also matches a primer that overhangs the end of the reference
//     (--min_overlap).
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
