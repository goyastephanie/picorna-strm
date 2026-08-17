//
// Two iterations of consensus assembly: reads are aligned to the selected
// reference, a rough consensus is called, reads are re-aligned to that
// consensus, and the final consensus is called from the second alignment.
//
// Optional per-alignment processing, driven by the sequencing mode:
//   - amplicon    : primer trimming (ivar trim)
//   - capture     : reserved; currently identical to metagenomic
//   - metagenomic : none
//

include { BWA_MEM_ALIGN as BWA_MEM_ALIGN_QUERY                                      } from '../modules/bwa_mem_align'
include { IVAR_CONSENSUS_BWA_ALIGN as IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY } from './ivar_consensus_bwa_align'
include { IVAR_CONSENSUS as BUILD_FINAL_CONSENSUS                                   } from '../modules/ivar_consensus'
include { MAKE_PRIMER_BED as PRIMER_BED_REF                                         } from '../modules/make_primer_bed'
include { MAKE_PRIMER_BED as PRIMER_BED_CON1                                        } from '../modules/make_primer_bed'
include { PRIMER_TRIM as PRIMER_TRIM_REF                                            } from '../modules/primer_trim'
include { PRIMER_TRIM as PRIMER_TRIM_CON1                                           } from '../modules/primer_trim'

// Alignments that failed the coverage/depth thresholds carry a sentinel BAM.
// ivar trim would choke on that placeholder, so it is routed around it and
// mixed back untouched, preserving the exact filename the downstream iVar step
// checks for.
def FAILED_BAM = 'FAILED.sorted.bam'

workflow CONSENSUS_ASSEMBLY {
    take:
    ch_reads         // channel: [ val(meta), path(reads) ]
    ch_ref           // channel: [ val(meta), val(ref_info), path(ref) ]
    use_mem2         // val: boolean
    do_primer_trim   // val: boolean
    primer_fwd       // val: forward primer sequence

    main:

    BWA_MEM_ALIGN_QUERY(
        ch_reads,
        ch_ref,
        use_mem2
    )

    //
    // ---- pass 1: alignment against the selected reference ----
    //
    ch_bam1 = BWA_MEM_ALIGN_QUERY.out.bam
    ch_ref1 = BWA_MEM_ALIGN_QUERY.out.ref

    if (do_primer_trim) {
        ch_bam1
            .join(ch_ref1, by: [0, 1])
            .branch { meta, ref_info, bam, bai, ref ->
                failed: bam.name == FAILED_BAM
                pass:   true
            }
            .set { ch_pt1 }

        ch_pt1.pass
            .multiMap { meta, ref_info, bam, bai, ref ->
                bam: [ meta, ref_info, bam, bai ]
                ref: [ meta, ref_info, ref ]
            }
            .set { ch_pt1_in }

        // locate the primer (python image), then trim by coordinate (iVar image)
        PRIMER_BED_REF ( ch_pt1_in.ref, primer_fwd )

        ch_pt1_in.bam
            .join(PRIMER_BED_REF.out.bed, by: [0, 1])
            .multiMap { meta, ref_info, bam, bai, bed ->
                bam: [ meta, ref_info, bam, bai ]
                bed: [ meta, ref_info, bed ]
            }
            .set { ch_pt1_trim }

        PRIMER_TRIM_REF ( ch_pt1_trim.bam, ch_pt1_trim.bed )

        ch_bam1 = PRIMER_TRIM_REF.out.bam
            .mix( ch_pt1.failed.map { meta, ref_info, bam, bai, ref -> [ meta, ref_info, bam, bai ] } )
    }


    // Re-synchronise: the optional steps above run asynchronously, so the BAM
    // order need not match the reference channel any more. iVar pairs its two
    // inputs positionally, so re-join by [meta, ref_info] and split again.
    ch_bam1
        .join(ch_ref1, by: [0, 1])
        .multiMap { meta, ref_info, bam, bai, ref ->
            bam: [ meta, ref_info, bam, bai ]
            ref: [ meta, ref_info, ref ]
        }
        .set { ch_pass1 }

    IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY(
        ch_pass1.bam,
        ch_pass1.ref,
        BWA_MEM_ALIGN_QUERY.out.reads,
        use_mem2
    )

    //
    // ---- pass 2: alignment against the first consensus ----
    //
    ch_bam2 = IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.bam
    ch_ref2 = IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.consensus

    if (do_primer_trim) {
        ch_bam2
            .join(ch_ref2, by: [0, 1])
            .branch { meta, ref_info, bam, bai, ref ->
                failed: bam.name == FAILED_BAM
                pass:   true
            }
            .set { ch_pt2 }

        ch_pt2.pass
            .multiMap { meta, ref_info, bam, bai, ref ->
                bam: [ meta, ref_info, bam, bai ]
                ref: [ meta, ref_info, ref ]
            }
            .set { ch_pt2_in }

        // locate the primer (python image), then trim by coordinate (iVar image)
        PRIMER_BED_CON1 ( ch_pt2_in.ref, primer_fwd )

        ch_pt2_in.bam
            .join(PRIMER_BED_CON1.out.bed, by: [0, 1])
            .multiMap { meta, ref_info, bam, bai, bed ->
                bam: [ meta, ref_info, bam, bai ]
                bed: [ meta, ref_info, bed ]
            }
            .set { ch_pt2_trim }

        PRIMER_TRIM_CON1 ( ch_pt2_trim.bam, ch_pt2_trim.bed )

        ch_bam2 = PRIMER_TRIM_CON1.out.bam
            .mix( ch_pt2.failed.map { meta, ref_info, bam, bai, ref -> [ meta, ref_info, bam, bai ] } )
    }


    ch_bam2
        .join(ch_ref2, by: [0, 1])
        .multiMap { meta, ref_info, bam, bai, ref ->
            bam: [ meta, ref_info, bam, bai ]
            ref: [ meta, ref_info, ref ]
        }
        .set { ch_pass2 }

    BUILD_FINAL_CONSENSUS (
        ch_pass2.bam,
        ch_pass2.ref
    )

    emit:
    final_consensus     = BUILD_FINAL_CONSENSUS.out.consensus
    initial_consensus   = IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.consensus
    bam                 = IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.bam   // channel: [ val(meta), val(ref_info), path(bam), path(bai) ]
    reads               = IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.reads // channel: [ val(meta), path(reads) ]
    init_covstats       = BWA_MEM_ALIGN_QUERY.out.covstats
    final_covstats      = IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.covstats
    // the exact alignment + reference the final consensus was called from,
    // ready for variant calling
    final_bam           = ch_pass2.bam
    final_ref           = ch_pass2.ref
}
