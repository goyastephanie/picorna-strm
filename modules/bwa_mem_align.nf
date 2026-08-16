process BWA_MEM_ALIGN {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_medium'
    container 'quay.io/epil02/revica-strm:0.0.5'

    input:
    tuple val(meta), path(fastq)
    tuple val(meta), val(ref_info), path(ref)
    val use_mem2

    output:
    tuple val(meta), val(ref_info), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    tuple val(meta), val(ref_info), path(ref),                                      emit: ref
    tuple val(meta), path(fastq),                                                   emit: reads
    tuple val(meta), path("*_failed_assembly.tsv"), optional: true,                 emit: failed_assembly
    tuple val(meta), val(ref_info), path("*_covstats.tsv"), optional: true,         emit: covstats

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ''
    def input = meta.single_end ? "${fastq}" : "${fastq[0]} ${fastq[1]}"
    def iter = task.ext.iter

    def min_coverage = task.ext.min_coverage
    def min_depth = task.ext.min_depth
    def bwa = use_mem2 ? "bwa-mem2" : "bwa"

    // this currently only supports ref paths to a fasta file, not a directory or subdirectory
    """
    # this flag combines 0x4 and 0x800, which means:
    # 1. read is unmapped
    # 2. read is supplementary alignment (chimeric, not representative alignment)
    # reads with either of these flags will be removed.
    FLAG=2052

    $bwa index $ref

    ## run bwa
    $bwa mem \
        $ref \
        $input \
        -t $task.cpus \
        | samtools view -bS -F \$FLAG -@ $task.cpus > ${prefix}

    samtools sort ${prefix} -@ ${task.cpus} -m 4G -o ${prefix}.bam
    samtools coverage -d 0 ${prefix}.bam > ${prefix}_depth.tsv

    mapped=\$(samtools view -c -F \$FLAG -@ $task.cpus ${prefix}.bam)

    prep_pandepth_output.py ${prefix}_depth.tsv $ref ${prefix}_covstats.tsv \\
        --extra_cols "reads_mapped_${iter}:\$mapped"

    coverage=\$(awk 'BEGIN {FS="\t"} NR>1 {print \$5}' "${prefix}_covstats.tsv")
    mean_depth=\$(awk 'BEGIN {FS="\t"} NR>1 {print \$4}' "${prefix}_covstats.tsv")

    # Check if thresholds are met
    if [ "\$(echo "\$coverage >= $min_coverage" | bc)" -eq 1 ] && [ "\$(echo "\$mean_depth >= $min_depth" | bc)" -eq 1 ]; then
        echo "alignment thresholds met!"
        samtools sort -@ ${task.cpus} -m 4G -o "${prefix}.sorted.bam" "${prefix}.bam"
        samtools index -@ ${task.cpus} "${prefix}.sorted.bam"
    else
        echo "alignment failed to meet thresholds! Depth: \$mean_depth, Coverage: \$coverage"
        rm *.bam
        touch FAILED.sorted.bam FAILED.sorted.bam.bai
        mv ${prefix}_covstats.tsv ${prefix}_failed_assembly.tsv
    fi
    """
}
