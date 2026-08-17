process BWA_MEM_ALIGN_DB {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/epil02/revica-strm:0.0.5'


    input:
    tuple val(meta), path(fastq)
    tuple path(db), path(db_indexed)
    val use_mem2

    output:
    tuple val(meta), path("*covstats.tsv"), emit: covstats
    tuple val(meta), path("*.bam") // emiting so we can save to output folder

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "$fastq" : "${fastq[0]} ${fastq[1]}"
    def bwa = use_mem2 ? "bwa-mem2" : "bwa"

    """
    # this flag combines 0x4 and 0x800, which means:
    # 1. read is unmapped
    # 2. read is supplementary alignment (chimeric, not representative alignment)
    # reads with either of these flags will be removed.
    FLAG=2052

    ## align and sort in a single pipe: writing an intermediate compressed BAM
    ## only to read it straight back costs a full zlib round trip plus two
    ## passes over hundreds of MB of the shared filesystem. `view -u` keeps the
    ## stream uncompressed on the way to sort. The sorted BAM is identical.
    $bwa mem \
        $db \
        $input \
        -t ${task.cpus} \
        | samtools view -u -F \$FLAG -@ ${task.cpus} \
        | samtools sort -@ ${task.cpus} -m 4G -o ${prefix}.bam

    samtools coverage -d 0 ${prefix}.bam > ${prefix}_depth.tsv

    ## replace abbreviated ref names in pandepth with originals from db
    ## we do this because BWA_MEM only records the alignment ref before the first space, which is
    ## usually just the acc number. We need the rest of the fasta header for downstream analyses

    prep_pandepth_output.py ${prefix}_depth.tsv $db ${prefix}_covstats.tsv
    """
}
