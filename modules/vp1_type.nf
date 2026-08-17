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
// Non-rhinovirus consensuses (enterovirus, coxsackievirus, echovirus,
// poliovirus) have no VP1 entry in this database at all and are correctly
// reported as "no_hit"; VP1 typing never gates assembly.
//
// -------------------------------------------------------------------------
// Why the hits are aggregated instead of taken straight from blast
// -------------------------------------------------------------------------
// blastn is a LOCAL aligner: it returns HSPs (high-scoring segment pairs),
// and one query/subject pair can produce several of them. An N counts as a
// mismatch during extension, so a run of Ns long enough to trip blast's
// X-drop cuts the alignment in two and blast re-seeds on the other side.
// Measured on this database with a real VP1 (OP342726.1_RV-A49, 849 nt),
// N runs spaced every 200 bases:
//
//     N run   total %N   reported HSP        with -max_hsps 1
//      <= 25    <= 9%    spans all 849 nt    typed correctly
//      >= 30      11%    breaks into ~200 nt no_hit (below --vp1_min_aln_len)
//
// So the old `-max_hsps 1 -perc_identity 72` combination silently lost VP1
// calls for consensuses with only ~10% Ns, purely because those Ns were
// clustered rather than scattered. Two things fix it, and they are the same
// two used by modules/check_contamination.nf:
//
//   * no -perc_identity filter, because blast's own identity counts every N
//     as a mismatch and would discard the hit before it is ever examined;
//
//   * all non-overlapping HSPs of a subject are summed, and identity is
//     computed ONLY over positions where both sequences call an unambiguous
//     A/C/G/T. Ns therefore neither create nor destroy identity.
//
// As a result `vp1_aln_len` is the number of unambiguous comparable positions
// (not blast's raw alignment length) and `vp1_identity` is the identity over
// exactly those positions. Both are stricter definitions than before, so the
// 72% / 300 nt defaults remain meaningful.
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
    def margin = params.vp1_ambiguous_margin
    """
    # Build the BLAST index inside the task so a stale/absent index in the
    # repository can never be used by accident.
    makeblastdb -in ${vp1_db} -dbtype nucl -out vp1db > /dev/null

    # -max_target_seqs is deliberately larger than the database: a low value is
    # applied early in the search and can drop the true best subject.
    # No -perc_identity (see header comment).
    blastn \\
        -query ${consensus} \\
        -db vp1db \\
        -max_target_seqs 500 \\
        -max_hsps 20 \\
        -evalue 1e-5 \\
        -outfmt "6 sseqid qstart qend sstart send bitscore qseq sseq" \\
        ${args} \\
        > raw_hits.tsv || true

    printf "sample\\tref_tag\\tref_acc\\tvp1_genotype\\tvp1_identity\\tvp1_aln_len\\tvp1_bitscore\\tvp1_subject\\tvp1_note\\n" \\
        > ${prefix}_vp1_type.tsv

    awk -F'\\t' \\
        -v minid="${min_id}" \\
        -v minaln="${min_aln}" \\
        -v margin="${margin}" \\
        -v sample="${meta.id}" \\
        -v reftag="${ref_info.tag}" \\
        -v refacc="${ref_info.acc}" '
    {
        s = \$1
        qst = \$2 + 0; qen = \$3 + 0; sst = \$4 + 0; sen = \$5 + 0
        if (qst > qen) { t = qst; qst = qen; qen = t }
        if (sst > sen) { t = sst; sst = sen; sen = t }

        # Keep an HSP only if it covers new territory in BOTH sequences. blast
        # lists the HSPs of a subject best-first, so the strongest alignment of
        # a region wins and a weaker duplicate of it is discarded.
        ov = 0
        for (j = 1; j <= nh[s]; j++) {
            if (qst <= qe[s, j] && qen >= qb[s, j]) { ov = 1; break }
            if (sst <= se[s, j] && sen >= sb[s, j]) { ov = 1; break }
        }
        if (ov) next
        j = ++nh[s]
        qb[s, j] = qst; qe[s, j] = qen
        sb[s, j] = sst; se[s, j] = sen
        bits[s] += \$6

        qseq = toupper(\$7); sseq = toupper(\$8)
        n = length(qseq)
        for (i = 1; i <= n; i++) {
            a = substr(qseq, i, 1); b = substr(sseq, i, 1)
            if (a != "A" && a != "C" && a != "G" && a != "T") continue
            if (b != "A" && b != "C" && b != "G" && b != "T") continue
            comp[s]++
            if (a == b) ident[s]++
        }
    }
    END {
        # best subject among those passing both thresholds, by summed bitscore
        best = ""; bestbits = -1
        for (s in comp) {
            if (comp[s] + 0 < minaln) continue
            if (100.0 * ident[s] / comp[s] < minid) continue
            if (bits[s] > bestbits) { bestbits = bits[s]; best = s }
        }

        if (best == "") {
            printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n", \\
                sample, reftag, refacc, "no_hit", "NA", "NA", "NA", "NA", \\
                "no_vp1_hit_above_threshold"
            exit
        }

        nf = split(best, ba, "_"); gt = ba[nf]
        pid = 100.0 * ident[best] / comp[best]

        # A close runner-up of a DIFFERENT genotype makes the call ambiguous.
        note = "ok"
        sec = ""; secbits = -1
        for (s in comp) {
            if (s == best) continue
            if (comp[s] + 0 < minaln) continue
            p = 100.0 * ident[s] / comp[s]
            if (p < minid) continue
            nf2 = split(s, sa, "_")
            if (sa[nf2] == gt) continue
            if (bits[s] > secbits) { secbits = bits[s]; sec = s }
        }
        if (sec != "") {
            sp = 100.0 * ident[sec] / comp[sec]
            if (pid - sp < margin) {
                nf3 = split(sec, sa2, "_")
                note = sprintf("ambiguous:%s(%.2f%%)", sa2[nf3], sp)
            }
        }

        printf "%s\\t%s\\t%s\\t%s\\t%.2f\\t%d\\t%.1f\\t%s\\t%s\\n", \\
            sample, reftag, refacc, gt, pid, comp[best], bits[best], best, note
    }' raw_hits.tsv >> ${prefix}_vp1_type.tsv
    """
}
