process BWA_MEM_INDEX {
    tag "${fasta}"
    label 'process_single'
    container 'quay.io/epil02/revica-strm:0.0.5'

    input: 
    path fasta

    output:
    tuple path(fasta), path("${fasta}*"), emit: indexed_fasta

    script:
    """
    bwa index $fasta
    """
}

    
