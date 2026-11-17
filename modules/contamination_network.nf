//
// Contamination network: which sample a flagged genome most plausibly came from.
//
// CHECK_CONTAMINATION says that two consensuses are near-identical; identity is
// symmetric, so it cannot say which of the two is the contaminant. This step
// reads that table together with the run summary and infers direction from
// abundance (bin/contamination_network.py):
//
//   1. possible_same_genome pairs are collapsed WITHIN each sample into one
//      observation (one virus assembled against two references). The
//      equivalence is never propagated to other samples.
//
//   2. every possible_contamination pair becomes an edge between two
//      collapsed sample-genome observations.
//
//   3. a low-abundance observation (<= --contam_contaminant_max_pct of the
//      sample's reads) linked to a high-abundance one (>= --contam_source_min_pct)
//      with a ratio of at least --contam_source_ratio is called
//      REMOVE_CONTAMINANT, and its partner KEEP_SOURCE. A second pass assigns an
//      indirect source when a connected component has exactly one
//      high-abundance candidate (source -> intermediate -> recipient).
//
// Everything that cannot be resolved safely -- two high-abundance samples
// sharing a genome, or a low-abundance genome with no source in the run -- is
// KEEP_UNRESOLVED, so a genuine co-infection is never removed on identity alone.
//
// Nothing is deleted: the outputs are decision tables. removal_manifest.tsv
// lists the sample/reference rows to drop, and cleaned_run_summary.tsv is the
// run summary without them.
//
process CONTAMINATION_NETWORK {
    tag "${params.run_name}"
    label 'process_low'
    container 'quay.io/biocontainers/pandas:2.2.1'

    input:
    path summary
    path comparisons

    output:
    path "contamination_results", emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    contamination_network.py \\
        --summary ${summary} \\
        --comparisons ${comparisons} \\
        --outdir contamination_results \\
        --metric ${params.contam_network_metric} \\
        --source-min-pct ${params.contam_source_min_pct} \\
        --contaminant-max-pct ${params.contam_contaminant_max_pct} \\
        --source-ratio ${params.contam_source_ratio} \\
        --high-conf-ratio ${params.contam_high_conf_ratio} \\
        --high-conf-identity ${params.contam_high_conf_identity}
    """
}
