//                                                                                 
// Consensus Assembly using BBMAP for alignment and iVar consensus for consensus calling  
//                                                              
                                                                                   
include { IVAR_CONSENSUS } from '../modules/ivar_consensus'                          
include { BWA_MEM_ALIGN } from '../modules/bwa_mem_align'
                                                                                   
workflow IVAR_CONSENSUS_BWA_ALIGN {                                                 
    take:                                                                          
    ch_bam                // channel: [ val(meta), val(ref_info), path(bam), path(bai) ]                             
    ch_ref                // channel: [ val(meta), val(ref_info), path(ref) ]                             
    ch_reads              // channel: [ val(meta), path(reads) ]                  
    use_mem2              // val:     use_mem2
                                                                                   
    main:                                                                          

    IVAR_CONSENSUS (
        ch_bam,
        ch_ref
    )

    ch_reads
        .join(IVAR_CONSENSUS.out.consensus)
        .multiMap { meta, reads, ref_info, consensus -> 
            reads:      [ meta, reads ]
            consensus:  [ meta, ref_info, consensus]
        }
        .set { ch_second_pass_input }

    
    BWA_MEM_ALIGN(
        ch_second_pass_input.reads,
        ch_second_pass_input.consensus,
        use_mem2
    )
    
    emit:
    bam         = BWA_MEM_ALIGN.out.bam       // channel: [ val(meta), val(ref_info), path(bam), path(bai) ]
    consensus   = BWA_MEM_ALIGN.out.ref       // channel: [ val(meta), val(ref_info), path(consensus) ]
    reads       = BWA_MEM_ALIGN.out.reads     // channel: [ val(meta), path(reads) ]
    covstats    = BWA_MEM_ALIGN.out.covstats
}
