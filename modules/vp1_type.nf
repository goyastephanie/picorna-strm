//
// Rhinovirus genotyping by VP1.
//
// Types a final consensus genome against a NUCLEOTIDE VP1 database
// (see bin/build_vp1_nt_db.py) using blastn. Rhinovirus types are defined by
// VP1 divergence, so this is the authoritative genotype call -- it can resolve
// genotype pairs that whole-genome mapping cannot separate (e.g. RV-A54 vs
// RV-A98), and it is reported independently of which genome reference the
// assembly happened to use.
//
// Non-rhinovirus consensuses (enterovirus, etc.) simply produce no hit above
// threshold and are reported as "no_hit".
//
process VP1_TYPE {
    tag "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    label 'process_single'
    container 'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0'

    input:
    tuple val(meta), val(ref_info), path(consensus)
    path vp1_db

    output:
    tuple val(meta), val(ref_info), path("*_vp1_type.tsv"), emit: genotype
    path "*_vp1_type.tsv",                                  emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${ref_info.acc}_${ref_info.tag}"
    def min_id = params.vp1_min_identity
    def min_aln = params.vp1_min_aln_len
    """
    # Build the BLAST index inside the task so a stale/absent index in the
    # repository can never be used by accident.
    makeblastdb -in ${vp1_db} -dbtype nucl -out vp1db > /dev/null

    blastn \\
        -query ${consensus} \\
        -db vp1db \\
        -max_target_seqs 5 \\
        -max_hsps 1 \\
        -perc_identity ${min_id} \\
        -outfmt "6 sseqid pident length bitscore" \\
        ${args} \\
        > raw_hits.tsv || true

    # Keep only hits with a long enough alignment, then take the best by bitscore.
    awk -F'\\t' -v m=${min_aln} '\$3 >= m' raw_hits.tsv | sort -t\$'\\t' -k4,4nr > hits.tsv

    printf "sample\\tref_tag\\tref_acc\\tvp1_genotype\\tvp1_identity\\tvp1_aln_len\\tvp1_bitscore\\tvp1_subject\\tvp1_note\\n" \\
        > ${prefix}_vp1_type.tsv

    if [ -s hits.tsv ]; then
        best=\$(head -1 hits.tsv)
        subj=\$(echo "\$best" | cut -f1)
        pid=\$(echo "\$best" | cut -f2)
        alen=\$(echo "\$best" | cut -f3)
        bits=\$(echo "\$best" | cut -f4)
        # genotype is the field after the last "_" of the subject id (ACC_RV-Axx)
        gt=\$(echo "\$subj" | awk -F'_' '{print \$NF}')

        # Flag a close runner-up of a DIFFERENT genotype: the call is then
        # ambiguous and deserves manual review.
        # NOTE: float comparison via awk, not bc -- the BLAST biocontainer is
        # minimal and has no bc.
        note="ok"
        second=\$(awk -F'\\t' -v g="\$gt" '{split(\$1,a,"_"); if (a[length(a)] != g) {print; exit}}' hits.tsv)
        if [ -n "\$second" ]; then
            spid=\$(echo "\$second" | cut -f2)
            sgt=\$(echo "\$second" | cut -f1 | awk -F'_' '{print \$NF}')
            close=\$(awk -v a="\$pid" -v b="\$spid" -v m=${params.vp1_ambiguous_margin} \\
                'BEGIN { print ((a - b) < m) ? 1 : 0 }')
            if [ "\$close" -eq 1 ]; then
                note="ambiguous:\${sgt}(\${spid}%)"
            fi
        fi

        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
            "${meta.id}" "${ref_info.tag}" "${ref_info.acc}" \\
            "\$gt" "\$pid" "\$alen" "\$bits" "\$subj" "\$note" \\
            >> ${prefix}_vp1_type.tsv
    else
        # No VP1 hit: either not a rhinovirus, or VP1 too incomplete/divergent.
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
            "${meta.id}" "${ref_info.tag}" "${ref_info.acc}" \\
            "no_hit" "NA" "NA" "NA" "NA" "no_vp1_hit_above_threshold" \\
            >> ${prefix}_vp1_type.tsv
    fi
    """
}
