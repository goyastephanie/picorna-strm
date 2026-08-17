process BWA_MEM_ALIGN {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_medium'
    container 'quay.io/epil02/revica-strm:0.0.5'

    input:
    tuple val(meta), path(fastq)
    tuple val(meta2), val(ref_info), path(ref)

    output:
    tuple val(meta), val(ref_info), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    tuple val(meta), val(ref_info), path(ref),                                      emit: ref
    // ref_info travels with the reads so that downstream re-joins can key on
    // [meta, ref_info]. A sample with a co-infection emits one item per selected
    // reference, all carrying the SAME meta, so meta alone is not a unique key
    // and `join` explicitly does not support duplicate keys.
    tuple val(meta), val(ref_info), path(fastq),                                    emit: reads
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

    // this currently only supports ref paths to a fasta file, not a directory or subdirectory
    """
    # this flag combines 0x4 and 0x800, which means:
    # 1. read is unmapped
    # 2. read is supplementary alignment (chimeric, not representative alignment)
    # reads with either of these flags will be removed.
    FLAG=2052

    bwa index $ref

    ## align and sort in a single pipe (see bwa_mem_align_db.nf): the previous
    ## version wrote an unsorted compressed BAM to disk, read it back, sorted it
    ## into a second BAM and later sorted THAT again into a third. The stream is
    ## kept uncompressed with `view -u` and sorted once; the result is identical.
    bwa mem \
        $ref \
        $input \
        -t $task.cpus \
        | samtools view -u -F \$FLAG -@ $task.cpus \
        | samtools sort -@ ${task.cpus} -m 4G -o ${prefix}.bam

    samtools coverage -d 0 ${prefix}.bam > ${prefix}_depth.tsv

    mapped=\$(samtools view -c -F \$FLAG -@ $task.cpus ${prefix}.bam)

    prep_pandepth_output.py ${prefix}_depth.tsv $ref ${prefix}_covstats.tsv \\
        --extra_cols "reads_mapped_${iter}:\$mapped"

    # covstats columns (see bin/prep_pandepth_output.py COV_COLS):
    # 1 #rname  2 startpos  3 endpos  4 numreads  5 covbases
    # 6 coverage(%)  7 meandepth  8 meanbaseq  9 meanmapq
    coverage=\$(awk 'BEGIN {FS="\t"} NR>1 {print \$6}' "${prefix}_covstats.tsv")
    mean_depth=\$(awk 'BEGIN {FS="\t"} NR>1 {print \$7}' "${prefix}_covstats.tsv")

    # Check if thresholds are met
    if [ "\$(echo "\$coverage >= $min_coverage" | bc)" -eq 1 ] && [ "\$(echo "\$mean_depth >= $min_depth" | bc)" -eq 1 ]; then
        echo "alignment thresholds met!"
        # already coordinate-sorted by the pipe above -- re-sorting it was a
        # no-op that cost a full pass and a third copy of the BAM
        mv "${prefix}.bam" "${prefix}.sorted.bam"
        samtools index -@ ${task.cpus} "${prefix}.sorted.bam"
    else
        echo "alignment failed to meet thresholds! Depth: \$mean_depth, Coverage: \$coverage"
        rm -f *.bam
        touch FAILED.sorted.bam FAILED.sorted.bam.bai
        mv ${prefix}_covstats.tsv ${prefix}_failed_assembly.tsv
    fi
    """
}
