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

import groovy.json.JsonOutput

// if INPUT not set
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// if db not set

if (!params.refs) { exit 1, "Reference database not specified!"}

// ---- sequencing mode -------------------------------------------------------
def VALID_MODES = ['amplicon', 'capture', 'metagenomic']
if (!(params.mode in VALID_MODES)) {
    exit 1, "Invalid --mode '${params.mode}'. Choose one of: ${VALID_MODES.join(', ')}"
}

// The mode sets the library-specific processing; the individual flags remain
// available as manual overrides.
def do_primer_trim = params.trim_primers || params.mode == 'amplicon'

if (do_primer_trim && !params.primer_fwd) {
    exit 1, "--mode amplicon (or --trim_primers) requires a forward primer (--primer_fwd)!"
}

log.info ""
log.info "  sequencing mode : ${params.mode}"
log.info "  primer trimming : ${do_primer_trim ? "yes (${params.primer_fwd})" : 'no'}"
log.info "  per-position VCF: ${params.call_variants ? 'yes' : 'no'}"
log.info "  iVar variants   : ${params.ivar_variants ? 'yes' : 'no'}"
log.info ""

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
include { IVAR_VARIANTS             } from './modules/ivar_variants'
include { BCFTOOLS_MPILEUP          } from './modules/bcftools_mpileup'
include { CHECK_CONTAMINATION       } from './modules/check_contamination'

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

//
// Software version registry.
//
// Every process runs in a pinned container, so the image tag IS the exact tool
// build that produced the results -- there is no way for the recorded version
// and the executed binary to disagree, which is not true of a `--version` call
// scraped at runtime. The registry therefore records, per process, the image
// Nextflow actually launched, plus the pipeline revision, the Nextflow build,
// the full command line, and an MD5 of each database so a re-run can be shown
// to have used the same sequences.
//
def md5_of(db_path) {
    try {
        def f = file(db_path)
        if (!f.exists()) return 'missing'
        if (f.size() > 500000000) return 'not_computed(file>500MB)'
        return java.security.MessageDigest.getInstance('MD5')
                   .digest(f.bytes).encodeHex().toString()
    } catch (Exception e) {
        return 'NA'
    }
}

def write_software_versions() {
    def y = new StringBuilder()
    y << "# picorna-strm software version registry\n"
    y << "# Container image tags are the authoritative record of tool versions.\n\n"
    y << "pipeline:\n"
    y << "  name: ${workflow.manifest.name ?: 'picorna-strm'}\n"
    y << "  version: ${workflow.manifest.version ?: 'unversioned'}\n"
    y << "  revision: ${workflow.revision ?: 'NA'}\n"
    y << "  commit: ${workflow.commitId ?: 'NA (run from a local directory)'}\n"
    y << "  projectDir: ${workflow.projectDir}\n\n"
    y << "run:\n"
    y << "  nextflow: ${workflow.nextflow.version}\n"
    y << "  nextflow_build: ${workflow.nextflow.build}\n"
    y << "  profile: ${workflow.profile}\n"
    y << "  container_engine: ${workflow.containerEngine ?: 'none'}\n"
    y << "  start: ${workflow.start}\n"
    y << "  complete: ${workflow.complete}\n"
    y << "  duration: ${workflow.duration}\n"
    y << "  success: ${workflow.success}\n"
    y << "  command_line: ${JsonOutput.toJson(workflow.commandLine)}\n\n"
    y << "databases:\n"
    y << "  refs: ${params.refs}\n"
    y << "  refs_md5: ${md5_of(params.refs)}\n"
    if (params.rhinovirus) {
        y << "  vp1_db: ${params.vp1_db}\n"
        y << "  vp1_db_md5: ${md5_of(params.vp1_db)}\n"
    }
    if (params.run_kraken2 && params.kraken2_db) {
        y << "  kraken2_db: ${params.kraken2_db}\n"
    }
    y << "\ncontainers:\n"
    def c = workflow.container
    if (c instanceof Map) {
        c.sort { it.key }.each { proc, img -> y << "  ${proc}: ${img ?: 'none'}\n" }
    } else {
        y << "  all_processes: ${c ?: 'none'}\n"
    }
    def dir = file("${params.output}/pipeline_info")
    dir.mkdirs()
    file("${params.output}/pipeline_info/software_versions.yml").text = y.toString()
}

workflow.onComplete {
    // never let bookkeeping mask the real outcome of the run
    try {
        def jsonStr = JsonOutput.toJson(params)
        file("${params.output}/params.json").text = JsonOutput.prettyPrint(jsonStr)
    } catch (Exception e) {
        log.warn "Could not write params.json: ${e.message}"
    }
    try {
        write_software_versions()
    } catch (Exception e) {
        log.warn "Could not write software_versions.yml: ${e.message}"
    }
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
                params.min_trimmed_reads,
                )

        ch_sample_input = FASTQ_TRIM_FASTP_MULTIQC.out.reads

        if (params.run_kraken2) {
            KRAKEN2 (
                    FASTQ_TRIM_FASTP_MULTIQC.out.reads,
                    file(params.kraken2_db),
                    params.kraken2_host_filter || params.save_kraken2_classified_reads || params.save_kraken2_unclassified_reads,
                    params.save_kraken2_classified_reads
                    )

                // actually swap in the host-depleted reads
                if (params.kraken2_host_filter) {
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
                params.use_mem2,
                do_primer_trim,
                params.primer_fwd
                )

        //
        // Variant calling against the sample's own final consensus, on the
        // same (primer-trimmed) alignment used to build it.
        //
        // Per-position pileup VCF against the sample's own consensus: the
        // record of intra-host / intra-genotype diversity.
        if (params.call_variants) {
            BCFTOOLS_MPILEUP (
                CONSENSUS_ASSEMBLY.out.final_bam,
                CONSENSUS_ASSEMBLY.out.final_ref
            )
        }

        // Optional filtered variant table with allele frequencies.
        if (params.ivar_variants) {
            IVAR_VARIANTS (
                CONSENSUS_ASSEMBLY.out.final_bam,
                CONSENSUS_ASSEMBLY.out.final_ref
            )
        }

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

        //
        // Cross-contamination: every final consensus of the run compared
        // against every other one. Reported as its own table rather than as a
        // column of the summary, because a summary row is written per assembly,
        // in parallel, long before the last consensus of the run exists --
        // contamination is a property of the pair, not of the sample.
        //
        if (params.check_contamination) {
            CHECK_CONTAMINATION (
                CONSENSUS_ASSEMBLY.out.final_consensus
                    .map { meta, ref_info, consensus -> consensus }
                    .collect()
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

