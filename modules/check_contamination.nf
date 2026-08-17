//
// Consensus-vs-consensus comparison of everything a run produced.
//
// It answers two different questions with the same measurement, because both
// have the same signature -- two consensuses that are near-identical over a
// long stretch -- and they are told apart only by whether the two came from the
// same sample:
//
//   * DIFFERENT samples -> `possible_contamination`. Index hopping, well-to-well
//     carry-over and a mis-pipetted library all produce two samples carrying the
//     same genome. Two genuinely separate infections, even within one outbreak,
//     almost always differ by at least a few bases.
//
//   * SAME sample -> `possible_same_genome`. Reference selection can pick two
//     references of two different genotypes for a single virus: if a sample
//     carries RV-A51, reads may also map well enough to an RV-A77 reference to
//     pass the coverage threshold, and the pipeline then builds two consensuses
//     that are really the same genome twice (VP1 typing usually gives away the
//     game by calling both A51). That is not a co-infection and not
//     contamination -- it is one virus masquerading as two, and one of the two
//     assemblies should be discarded. A true co-infection of two genotypes of
//     one species sits around 75-85% identity genome-wide, nowhere near the
//     threshold.
//
// Only consensuses of the SAME species are compared (--contam_intraspecies_only,
// on by default): the species is the tag with its trailing serotype number
// removed, so RV-A51 and RV-A77 are compared, RV-A51 and RV-C11 are not.
// Between species the identity is far below any threshold used here, so the
// restriction costs nothing and keeps the table free of 5'UTR-only hits.
// NOTE: for enteroviruses this prefix is the serotype label, not the ICTV
// species -- E11, CV-B3 and EV-B69 are all Enterovirus B but carry different
// prefixes, so they are not compared to each other. Set
// --contam_intraspecies_only false to compare everything.
//
// Three details decide whether the reported number means anything:
//
//   * consensuses have different lengths, so only the aligned overlap is
//     compared and its size is reported (`comparable_sites`). 100% identity
//     over 300 nt is not evidence of anything; over 5000 nt it is.
//
//   * Ns and IUPAC ambiguity codes are excluded from BOTH the numerator and the
//     denominator. A position that is N in one genome says nothing about
//     whether the two share it; counting it as a difference would hide the
//     signal in exactly the low-coverage assemblies where it matters most, and
//     counting it as a match would invent identity. This also makes the
//     same-sample check work as intended, since two assemblies of one virus
//     differ mostly by which regions each of them left ambiguous.
//
//   * blastn is run WITHOUT -perc_identity. That filter uses blast's own
//     identity, in which every N is a mismatch, so a genuine contaminant with
//     20-30% Ns would be dropped before it is ever examined. Non-overlapping
//     HSPs are summed instead, which recovers the segments between long N runs.
//
// A flag is a lead, not a verdict. For `possible_contamination`, confirm by
// looking for the partner's alleles as low-frequency variants in the VCF of the
// suspected recipient. For `possible_same_genome`, compare the two VP1 calls and
// the coverage of each assembly in the summary.
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
    def same_genome_identity = params.contam_same_genome_identity
    def intraspecies_only = params.contam_intraspecies_only ? 1 : 0
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

        awk -F'\\t' \\
            -v minc="${min_comparable}" \\
            -v flagid="${flag_identity}" \\
            -v sameid="${same_genome_identity}" \\
            -v intra="${intraspecies_only}" '
        {
            q = \$1; s = \$2
            if (q == s) next
            key = (q < s) ? q "\\001" s : s "\\001" q

            # accumulate HSPs from one of the two reciprocal blast directions only
            if (!(key in dir)) {
                split(q, qa, "|"); split(s, sa, "|")

                # species = tag with the trailing serotype number removed
                # (RV-A51 -> RV-A, EV-D68 -> EV-D, E11 -> E)
                spa = qa[3]
                if (match(spa, /[0-9]+[A-Za-z]*\$/)) spa = substr(spa, 1, RSTART - 1)
                if (spa == "") spa = qa[3]
                spb = sa[3]
                if (match(spb, /[0-9]+[A-Za-z]*\$/)) spb = substr(spb, 1, RSTART - 1)
                if (spb == "") spb = sa[3]
                if (intra == 1 && spa != spb) { dir[key] = "SKIP"; next }

                dir[key] = q
                same[key] = (qa[1] == sa[1]) ? 1 : 0
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
                if (same[key] == 1) {
                    # one sample, two assemblies: never contamination
                    flag = (pid >= sameid) ? "possible_same_genome" : "ok"
                    rank = (pid >= sameid) ? 2 : 3
                } else {
                    flag = (pid >= flagid) ? "possible_contamination" : "ok"
                    rank = (pid >= flagid) ? 1 : 3
                }
                printf "%d\\t%s\\t%s\\t%d\\t%d\\t%.4f\\t%s\\n", rank, la[key], lb[key], c, c - ident[key], pid, flag
            }
        }' raw_pairs.tsv \\
        | sort -t"\$(printf '\\t')" -k1,1n -k8,8nr -k2,2 \\
        | cut -f2- \\
        >> ${prefix}_contamination.tsv
    fi

    n_contam=\$(awk -F'\\t' 'NR>1 && \$8 == "possible_contamination"' ${prefix}_contamination.tsv | wc -l)
    n_same=\$(awk -F'\\t' 'NR>1 && \$8 == "possible_same_genome"' ${prefix}_contamination.tsv | wc -l)
    echo "[check_contamination] \$n_seq consensuses, \$n_contam possible contamination, \$n_same possible duplicate assembly"
    """
}
