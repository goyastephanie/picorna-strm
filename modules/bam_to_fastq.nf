process BAM_TO_FASTQ {

    container 'quay.io/epil02/revica-strm:0.0.6'
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_low'

    input:
    tuple val(meta), val(ref_info), path(bam), path(bai)
    val sra_proper_pair

    output:
    path "*_SRA.fastq.gz", optional: true, emit: fastq

    script:
    def prefix = task.ext.prefix ?: ''

    def output = ""
    def prep = ""
    def inbam = "$bam"
    def cleanup = ""

    // single end, dump all reads into a single file.
    if (meta.single_end) {
        output = "-0 ${prefix}_SRA.fastq.gz"
    }

    // paired end, dump all reads into a single file, regardless of proper pairing.
    if (!sra_proper_pair & !meta.single_end) {
        output = "-o ${prefix}_SRA.fastq.gz"
    }

    // paired end, dump only properly-paired reads into a single file.
    if (sra_proper_pair & !meta.single_end) {
        output = "-o ${prefix}_SRA.fastq.gz -s /dev/null"
        prep = "samtools sort -o ${prefix}_name_sorted.bam -m 4G -N -@ ${task.cpus} $bam" 
        inbam = "${prefix}_name_sorted.bam"
        cleanup = "rm ${prefix}_name_sorted.bam"
    }

    """

    if [[ \$(basename "$bam") = "FAILED.sorted.bam" ]]; then
        echo "Skipping bam to fastq conversion; alignment with ${prefix} failed depth/coverage previously"
        exit 0 # shouldn't cause fail if the outputs are optional
    fi

    $prep

    samtools fastq \\
    $inbam \\
    -F 4 \\
    -N \\
    -@ ${task.cpus} \\
    $output

    $cleanup

    """
}
