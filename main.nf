#!/usr/bin/env nextflow
/*
   ========================================================================================
   picorna-strm
   ========================================================================================
   Based on revica-strm (https://github.com/greninger-lab/revica-strm),
   adapted for picornavirus genomes.

Author:
Eli Piliper <epil02@uw.edu>
Jaydee Sereewit <aseree@uw.edu>
Stephanie Goya <sgoya@uw.edu>
Alex L Greninger <agrening@uw.edu>
UW Medicine | Virology
Department of Laboratory Medicine and Pathology
University of Washington
LICENSE: GNU
----------------------------------------------------------------------------------------
 */

// if INPUT not set
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// if db not set

if (!params.refs) { exit 1, "Reference database not specified!"}

// optional rhinovirus VP1 genotyping module
if (params.rhinovirus) {
    if (!params.vp1_db) { exit 1, "--rhinovirus requires a VP1 nucleotide database (--vp1_db)!" }
    if (!file(params.vp1_db).exists()) { exit 1, "VP1 database not found: ${params.vp1_db}" }
}
//
// SUBWORKFLOWS
//
include { INPUT_CHECK               } from './subworkflows/input_check'
include { REFERENCE_PREP            } from './subworkflows/reference_prep'
include { CONSENSUS_ASSEMBLY        } from './subworkflows/consensus_assembly'
include { FASTQ_TRIM_FASTP_MULTIQC  } from './subworkflows/fastq_trim_fastp_multiqc'

//
// MODULES
//
include { SEQTK_SAMPLE              } from './modules/seqtk_sample'
include { SUMMARY                   } from './modules/summary'
include { KRAKEN2                   } from './modules/kraken2'
include { BAM_TO_FASTQ              } from './modules/bam_to_fastq'
include { VP1_TYPE                  } from './modules/vp1_type'

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
log.info "                                                                                    "
log.info " ██████╗ ██╗ ██████╗ ██████╗ ██████╗ ███╗   ██╗ █████╗       ███████╗████████╗██████╗ ███╗   ███╗ "
log.info " ██╔══██╗██║██╔════╝██╔═══██╗██╔══██╗████╗  ██║██╔══██╗      ██╔════╝╚══██╔══╝██╔══██╗████╗ ████║ "
log.info " ██████╔╝██║██║     ██║   ██║██████╔╝██╔██╗ ██║███████║█████╗███████╗   ██║   ██████╔╝██╔████╔██║ "
log.info " ██╔═══╝ ██║██║     ██║   ██║██╔══██╗██║╚██╗██║██╔══██║╚════╝╚════██║   ██║   ██╔══██╗██║╚██╔╝██║ "
log.info " ██║     ██║╚██████╗╚██████╔╝██║  ██║██║ ╚████║██║  ██║      ███████║   ██║   ██║  ██║██║ ╚═╝ ██║ "
log.info " ╚═╝     ╚═╝ ╚═════╝ ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝      ╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝     ╚═╝ "
log.info "        reference-based assembly + VP1 genotyping for picornaviruses                 "
log.info "                                                                                    "

// Save run params to file
import groovy.json.JsonOutput

workflow.onComplete {
    jsonStr = JsonOutput.toJson(params)
    file("${params.output}/params.json").text = JsonOutput.prettyPrint(jsonStr)
}

workflow {

        INPUT_CHECK (
                ch_input
                )

        FASTQ_TRIM_FASTP_MULTIQC (
                INPUT_CHECK.out.reads,
                params.adapter_trimming,
                params.save_trimmed_fail,
                params.save_merged,
                params.skip_fastp,
                params.min_trimmed_reads,
                )

        ch_sample_input = FASTQ_TRIM_FASTP_MULTIQC.out.reads

        if (params.run_kraken2) {
            KRAKEN2 (
                    FASTQ_TRIM_FASTP_MULTIQC.out.reads,
                    file(params.kraken2_db),
                    params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter,
                    params.kraken2_variants_host_filter || params.kraken2_assembly_host_filter
                    )

                if (params.kraken2_variants_host_filter) {
                    ch_sample_input = KRAKEN2.out.unclassified_reads_fastq
                }
        }

    if (params.sample) {
        SEQTK_SAMPLE (
                ch_sample_input,
                params.sample
                )
            ch_ref_prep_input = SEQTK_SAMPLE.out.reads
    } else {
        ch_ref_prep_input = ch_sample_input
    } 

    if (!params.skip_consensus) { 
        REFERENCE_PREP (
                ch_ref_prep_input,
                file(params.refs),
                params.use_mem2
                ) 

        CONSENSUS_ASSEMBLY (
                REFERENCE_PREP.out.reads,
                REFERENCE_PREP.out.ref,
                params.use_mem2
                )

        BAM_TO_FASTQ(
            CONSENSUS_ASSEMBLY.out.bam,
            params.sra_proper_pair
        )

        //
        // Optional: rhinovirus genotyping by VP1.
        //
        // Deliberately decoupled from reference selection: the assembly still
        // uses the best whole-genome reference, and VP1 provides the
        // authoritative genotype call on top of the finished consensus.
        //
        if (params.rhinovirus) {
            VP1_TYPE (
                CONSENSUS_ASSEMBLY.out.final_consensus,
                file(params.vp1_db)
            )

            VP1_TYPE.out.tsv
                .collectFile(
                    storeDir: "${params.output}",
                    name: "${params.run_name}_vp1_genotypes.tsv",
                    keepHeader: true,
                    sort: true
                )
        }

        ch_summary_in = FASTQ_TRIM_FASTP_MULTIQC.out.trim_log
            .combine(
                CONSENSUS_ASSEMBLY.out.final_consensus
                    .join(CONSENSUS_ASSEMBLY.out.init_covstats, by: [0,1])
                    .join(CONSENSUS_ASSEMBLY.out.final_covstats, by: [0,1])
                , by: 0
            ).map {meta, trim_log, ref_info, consensus, init_covstats, final_covstats -> [meta, ref_info, trim_log, consensus, init_covstats, final_covstats]}

        SUMMARY (
                ch_summary_in
                )

            SUMMARY.out.summary
            .collectFile(storeDir: "${params.output}", name:"${params.run_name}_summary.tsv", keepHeader: true, sort: true)

            SUMMARY.out.ready_to_concat
            .collect()
            .map { it -> true }
        .set { all_summaries_done }

    }
}

