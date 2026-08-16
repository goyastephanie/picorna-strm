//
// 2 iterations of consensus assembly using BWA MEM for alignment and iVar consensus for consensus calling
//

include { BWA_MEM_ALIGN as BWA_MEM_ALIGN_QUERY                                     } from '../modules/bwa_mem_align'
include { IVAR_CONSENSUS_BWA_ALIGN as IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY } from './ivar_consensus_bwa_align'
include { IVAR_CONSENSUS as BUILD_FINAL_CONSENSUS                                    } from '../modules/ivar_consensus'
include { PICARD_MARKDUPLICATES as DEDUP_REF                                         } from '../modules/picard_markduplicates'
include { PICARD_MARKDUPLICATES as DEDUP_CON1                                        } from '../modules/picard_markduplicates'

workflow CONSENSUS_ASSEMBLY {
    take:
    ch_reads    // channel: [ val(meta), path(reads) ]
    ch_ref      // channel: [ val(meta), val(ref_info), path(ref) ]
    use_mem2    // val use_mem2

    main:

    BWA_MEM_ALIGN_QUERY(
        ch_reads,
        ch_ref,
        use_mem2
    )

    //
    // Optionally remove PCR duplicates from the reference alignment before
    // building the initial consensus. Recommended for hybridization-capture
    // libraries; must stay OFF for amplicon data (reads share start
    // coordinates by design and would be wrongly collapsed).
    //
    // Notes:
    //   - Failed alignments carry a sentinel "FAILED.sorted.bam"; Picard would
    //     choke on that dummy file, so we branch those out, dedup only real
    //     BAMs, then mix the sentinels back untouched.
    //   - Picard runs asynchronously, so its output order need not match the
    //     reference channel. iVar pairs its two inputs positionally, so we
    //     re-join by [meta, ref_info] and split back with multiMap to keep the
    //     BAM and reference channels in lockstep.
    //
    ch_init_bam = BWA_MEM_ALIGN_QUERY.out.bam
    ch_init_ref = BWA_MEM_ALIGN_QUERY.out.ref
    if (params.remove_duplicates) {
        BWA_MEM_ALIGN_QUERY.out.bam
            .branch { meta, ref_info, bam, bai ->
                failed: bam.name == 'FAILED.sorted.bam'
                pass:   true
            }
            .set { ch_ref_bam_split }

        DEDUP_REF ( ch_ref_bam_split.pass )

        DEDUP_REF.out.bam
            .mix( ch_ref_bam_split.failed )
            .join( BWA_MEM_ALIGN_QUERY.out.ref, by: [0, 1] )
            .multiMap { meta, ref_info, bam, bai, ref ->
                bam: [ meta, ref_info, bam, bai ]
                ref: [ meta, ref_info, ref ]
            }
            .set { ch_init }

        ch_init_bam = ch_init.bam
        ch_init_ref = ch_init.ref
    }

    IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY(
        ch_init_bam,
        ch_init_ref,
        BWA_MEM_ALIGN_QUERY.out.reads,
        use_mem2
    )

    //
    // Optionally remove PCR duplicates from the alignment to the initial
    // consensus before building the final (delivered) consensus.
    //
    ch_final_bam       = IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.bam
    ch_final_consensus = IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.consensus
    if (params.remove_duplicates) {
        IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.bam
            .branch { meta, ref_info, bam, bai ->
                failed: bam.name == 'FAILED.sorted.bam'
                pass:   true
            }
            .set { ch_con1_bam_split }

        DEDUP_CON1 ( ch_con1_bam_split.pass )

        DEDUP_CON1.out.bam
            .mix( ch_con1_bam_split.failed )
            .join( IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.consensus, by: [0, 1] )
            .multiMap { meta, ref_info, bam, bai, consensus ->
                bam:       [ meta, ref_info, bam, bai ]
                consensus: [ meta, ref_info, consensus ]
            }
            .set { ch_final }

        ch_final_bam       = ch_final.bam
        ch_final_consensus = ch_final.consensus
    }

    BUILD_FINAL_CONSENSUS (
        ch_final_bam,
        ch_final_consensus
    )

    emit:
    final_consensus     = BUILD_FINAL_CONSENSUS.out.consensus
    initial_consensus   = IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.consensus
    bam                 = IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.bam       // channel: [ val(meta), val(ref_info), path(bam), path(bai) ]
    reads               = IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.reads     // channel: [ val(meta), path(reads) ]
    init_covstats       = BWA_MEM_ALIGN_QUERY.out.covstats
    final_covstats      = IVAR_CONSENSUS_BWA_MEM_ALIGN_INITIAL_ASSEMBLY.out.covstats
}
