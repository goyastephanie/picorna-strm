process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0'
    tag "OUTDIR: ${output_dir}"
    label = 'process_high' 
    errorStrategy = 'ignore'

    input:
    val run_name
    val output_dir 
    path(fastp_files)

    output:
    path "${run_name}_multiqc.html", emit: report

    """
    echo "FastP files:"
    echo "${fastp_files}"
    
    multiqc -f \\
        -o . \\
        -n ${run_name}_multiqc.html \\
        -m fastp \\
        ${fastp_files}
    """
}
