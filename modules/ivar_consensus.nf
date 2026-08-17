process IVAR_CONSENSUS {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    // samtools mpileup | ivar consensus are BOTH single-threaded (this module
    // never references task.cpus). 'process_high' would reserve every core of the
    // machine to use one, serialising the whole consensus stage.
    label 'process_single'
    container 'quay.io/biocontainers/ivar:1.4.4--h077b44d_0'

    input:
    tuple val(meta), val(ref_info), path(bam), path(bai)
    tuple val(meta2), val(ref_info2), path(ref)

    output:
    tuple val(meta), val(ref_info), path("*.fa"),       optional: true, emit: consensus
    tuple val(meta), val(ref_info), path("*.qual.txt"), optional: true, emit: qual
    tuple val(meta), val(ref_info), path("*.mpileup"),  optional: true, emit: mpileup

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: ''
    """

    if [[ \$(basename "$bam") = "FAILED.sorted.bam" ]]; then
        echo "Skipping $prefix consensus with $ref; failed depth/coverage previously"
        rm *.fa
        exit 0 # shouldn't cause fail if the outputs are optional
    fi

    samtools \\
        mpileup \\
        --reference $ref \\
        $args2 \\
        $bam \\
        | ivar \\
            consensus \\
            $args \\
            -p $prefix \\
            -i $prefix

    # get rid of linebreaks except the header line
    awk '/^>/ {printf "%s\\n", \$0; next} {printf "%s", \$0} END {print ""}' ${prefix}.fa > ${prefix}_temp.fa
    # removing leading Ns
    sed '/^>/!s/^N\\+//' ${prefix}_temp.fa > ${prefix}_temp_frontNtrimmed.fa
    # remove trailing Ns
    sed '/^>/!s/N\\+\$//' ${prefix}_temp_frontNtrimmed.fa > ${prefix}.fa
    
    rm ${prefix}_temp.fa
    rm ${prefix}_temp_frontNtrimmed.fa
    """
}
