include { FASTP   } from '../modules/fastp'
include { MULTIQC } from '../modules/multiqc'

import groovy.json.JsonSlurper

def getFastpReadsAfterFiltering(json_file) {
    def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
    return json['after_filtering']['total_reads'].toInteger()
}

workflow FASTQ_TRIM_FASTP_MULTIQC {

    take:
    ch_reads                // channel: [ val(meta), path(reads)  ]
    ch_adapter_fasta        // channel: [ path(fasta) ]
    val_save_trimmed_fail   // value: boolean
    val_save_merged         // value: boolean
    val_min_read_trim_limit // value: int 

    main:
    // fastp always runs: every downstream step assumes quality-filtered reads,
    // and MULTIQC/SUMMARY both depend on its reports.
    FASTP (
        ch_reads,
        ch_adapter_fasta,
        val_save_trimmed_fail,
        val_save_merged
    )
    ch_trim_reads        = FASTP.out.reads
    ch_trim_json         = FASTP.out.json
    ch_trim_html         = FASTP.out.html
    ch_trim_log          = FASTP.out.log
    ch_trim_reads_fail   = FASTP.out.reads_fail
    ch_trim_reads_merged = FASTP.out.reads_merged

    // filter first, THEN map: a map whose closure has no else branch returns
    // null for the samples it rejects, and that null is emitted into the
    // channel and breaks the next process that destructures it
    ch_trim_reads
        .join(ch_trim_json)
        .filter { meta, reads, json ->
            def n = getFastpReadsAfterFiltering(json)
            // a dropped sample produces no assembly, no summary row and no entry
            // under fail/ -- say so on the console, or it just vanishes
            if (n < val_min_read_trim_limit) {
                log.warn "${meta.id}: ${n} reads after fastp, below --min_trimmed_reads ${val_min_read_trim_limit}; sample excluded from assembly"
            }
            return n >= val_min_read_trim_limit
        }
        .map { meta, reads, json -> [ meta, reads ] }
        .set { ch_trim_reads }

    MULTIQC (
        params.run_name,
        file("${params.output}").toAbsolutePath().toString(),
        ch_trim_html.map{meta, file -> file.parent.toAbsolutePath() }.collect(),
    )

    emit:
    reads             = ch_trim_reads         // channel: [ val(meta), path(reads) ]
    trim_json         = ch_trim_json          // channel: [ val(meta), path(json) ]
    trim_html         = ch_trim_html          // channel: [ val(meta), path(html) ]
    trim_log          = ch_trim_log           // channel: [ val(meta), path(log) ]
    trim_reads_fail   = ch_trim_reads_fail    // channel: [ val(meta), path(fastq.gz) ]
    trim_reads_merged = ch_trim_reads_merged  // channel: [ val(meta), path(fastq.gz) ]
    multiqc_report    = MULTIQC.out.report
}
