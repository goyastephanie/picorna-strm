process BWA_MEM_INDEX {
    tag "${fasta}"
    label 'process_single'
    container 'quay.io/epil02/revica-strm:0.0.5'

    input: 
    path fasta
    val use_mem2

    output:
    tuple path(fasta), path("${fasta}*"), emit: indexed_fasta

    script:
    def bwa = use_mem2 ? "bwa-mem2" : "bwa"

    """
    $bwa index $fasta
    """
}

    
