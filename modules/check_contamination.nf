//
// Cross-contamination check: all-vs-all comparison of the final consensuses
// produced by a run.
//
// Index hopping, well-to-well carry-over and a mis-pipetted library all produce
// the same signature: two genomes from DIFFERENT samples that are *identical*
// over a long stretch. Two genuinely separate infections, even in the same
// outbreak, almost always differ by at least a few bases.
//
// Three details decide whether the number reported here means anything:
//
//   * consensuses have different lengths, so only the aligned overlap is
//     compared and its size is reported (`comparable_sites`). 100% identity
//     over 300 nt is not evidence of anything; over 5000 nt it is.
//
//   * Ns and IUPAC ambiguity codes are excluded from BOTH the numerator and the
//     denominator. A position that is N in one genome says nothing about
//     whether the two samples share it; counting it as a difference would hide
//     contamination in exactly the low-coverage samples where contamination
//     matters most, and counting it as a match would invent identity.
//
//   * blastn is run WITHOUT -perc_identity. That filter uses blast's own
//     identity, in which every N is a mismatch, so a genuine contaminant with
//     20-30% Ns would be dropped before it is ever examined. Non-overlapping
//     HSPs are summed instead, which recovers the segments between long N runs.
//
// Pairs from the same sample (a real co-infection assembled against two
// references) are skipped by construction.
//
// A flag is a lead, not a verdict. Confirm it by looking for the partner's
// alleles as low-frequency variants in the VCF of the suspected recipient.
//
process CHECK_CONTAMINATION {
    tag "${params.run_name}"
    label 'process_low'
    container 'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0'

    input:
    path consensus

    output:
    path "*_contamination.tsv", emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${params.run_name}"
    def min_comparable = params.contam_min_comparable
    def flag_identity = params.contam_flag_identity
    """
    # One FASTA holding every consensus, with the sample / accession / tag
    # recovered from the file name and encoded in the header.
    for f in *_consensus_final.fa; do
        base=\$(basename "\$f" _consensus_final.fa)
        tag=\${base##*_}
        rest=\${base%_*}
        acc=\${rest##*_}
        sample=\${rest%_*}
        awk -v s="\$sample" -v a="\$acc" -v t="\$tag" \\
            '/^>/ { n++; print ">" s "|" a "|" t "|" n; next } { print }' "\$f"
    done > all_consensus.fa

    n_seq=\$(grep -c '^>' all_consensus.fa || true)

    printf "sample_a\\ttag_a\\tsample_b\\ttag_b\\tcomparable_sites\\tdifferences\\tpct_identity\\tflag\\n" \\
        > ${prefix}_contamination.tsv

    if [ "\$n_seq" -ge 2 ]; then
        makeblastdb -in all_consensus.fa -dbtype nucl -out ctdb > /dev/null

        # NOTE: no -perc_identity here on purpose (see header comment).
        blastn \\
            -query all_consensus.fa \\
            -db ctdb \\
            -outfmt "6 qseqid sseqid qstart qend sstart send qseq sseq" \\
            -max_target_seqs 5000 \\
            -max_hsps 20 \\
            -evalue 1e-20 \\
            > raw_pairs.tsv || true

        awk -F'\\t' -v minc="${min_comparable}" -v flagid="${flag_identity}" '
        {
            q = \$1; s = \$2
            if (q == s) next
            key = (q < s) ? q "\\001" s : s "\\001" q

            # accumulate HSPs from one of the two reciprocal blast directions only
            if (!(key in dir)) {
                split(q, qa, "|"); split(s, sa, "|")
                if (qa[1] == sa[1]) { dir[key] = "SKIP"; next }   # same sample
                dir[key] = q
                la[key] = qa[1] "\\t" qa[3]
                lb[key] = sa[1] "\\t" sa[3]
            }
            if (dir[key] == "SKIP" || dir[key] != q) next

            qst = \$3 + 0; qen = \$4 + 0; sst = \$5 + 0; sen = \$6 + 0
            if (qst > qen) { t = qst; qst = qen; qen = t }
            if (sst > sen) { t = sst; sst = sen; sen = t }

            # keep an HSP only if it covers new territory in both sequences;
            # blast returns HSPs best-first, so the strongest one wins a region
            ov = 0
            for (j = 1; j <= nhsp[key]; j++) {
                if (qst <= qe[key, j] && qen >= qb[key, j]) { ov = 1; break }
                if (sst <= se[key, j] && sen >= sb[key, j]) { ov = 1; break }
            }
            if (ov) next
            j = ++nhsp[key]
            qb[key, j] = qst; qe[key, j] = qen
            sb[key, j] = sst; se[key, j] = sen

            qseq = toupper(\$7); sseq = toupper(\$8)
            n = length(qseq)
            for (i = 1; i <= n; i++) {
                a = substr(qseq, i, 1); b = substr(sseq, i, 1)
                if (a != "A" && a != "C" && a != "G" && a != "T") continue
                if (b != "A" && b != "C" && b != "G" && b != "T") continue
                comp[key]++
                if (a == b) ident[key]++
            }
        }
        END {
            for (key in dir) {
                if (dir[key] == "SKIP") continue
                c = comp[key] + 0
                if (c < minc) continue
                pid = 100.0 * ident[key] / c
                flag = (pid >= flagid) ? "possible_contamination" : "ok"
                printf "%s\\t%s\\t%d\\t%d\\t%.4f\\t%s\\n", la[key], lb[key], c, c - ident[key], pid, flag
            }
        }' raw_pairs.tsv \\
        | sort -t"\$(printf '\\t')" -k7,7nr -k1,1 \\
        >> ${prefix}_contamination.tsv
    fi

    n_flag=\$(awk -F'\\t' 'NR>1 && \$8 == "possible_contamination"' ${prefix}_contamination.tsv | wc -l)
    echo "[check_contamination] \$n_seq consensuses, \$n_flag flagged pair(s)"
    """
}
