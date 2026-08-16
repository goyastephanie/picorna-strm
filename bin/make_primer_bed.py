#!/usr/bin/env python3
"""
Locate a primer sequence in a reference FASTA and emit a BED file for ivar trim.

Handles IUPAC ambiguity codes in the primer (e.g. Y = C/T) and tolerates a
configurable number of mismatches, because the primer-binding site of a given
sample is rarely identical to the reference.

Both strands are searched. If the primer is not found, an EMPTY BED is written
and the script still exits 0: ivar trim then simply has nothing to trim, which
is the correct behaviour for a reference whose 5' end is missing.
"""

import argparse
import sys

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
COMP = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn",
                     "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def find_matches(seq, primer, max_mm):
    """Return list of (start0, end, n_mismatch) for primer matches in seq."""
    seq = seq.upper()
    primer = primer.upper()
    plen = len(primer)
    allowed = [IUPAC.get(b, "ACGT") for b in primer]
    hits = []
    for i in range(0, len(seq) - plen + 1):
        mm = 0
        for j in range(plen):
            if seq[i + j] not in allowed[j]:
                mm += 1
                if mm > max_mm:
                    break
        if mm <= max_mm:
            hits.append((i, i + plen, mm))
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="reference FASTA")
    ap.add_argument("--primer", required=True, help="primer sequence (IUPAC allowed)")
    ap.add_argument("--out", required=True, help="output BED")
    ap.add_argument("--name", default="primer_fwd", help="primer name for the BED")
    ap.add_argument("--max_mismatch", type=int, default=4,
                    help="mismatches tolerated when locating the primer (default: 4)")
    ap.add_argument("--all_matches", action="store_true",
                    help="report every match instead of only the best one per strand")
    args = ap.parse_args()

    rows = []
    for name, seq in read_fasta(args.ref):
        for strand, s in (("+", seq), ("-", revcomp(seq))):
            hits = find_matches(s, args.primer, args.max_mismatch)
            if not hits:
                continue
            if not args.all_matches:
                hits = [sorted(hits, key=lambda h: (h[2], h[0]))[0]]
            for start, end, mm in hits:
                if strand == "-":
                    # translate coordinates back onto the forward strand
                    start, end = len(seq) - end, len(seq) - start
                rows.append((name, start, end, f"{args.name}_{strand}", mm, strand))

    with open(args.out, "w") as out:
        for name, start, end, pname, mm, strand in rows:
            out.write(f"{name}\t{start}\t{end}\t{pname}\t{60 - mm}\t{strand}\n")

    if rows:
        for r in rows:
            sys.stderr.write(
                f"[make_primer_bed] {r[0]}:{r[1]}-{r[2]} ({r[5]}) "
                f"{r[3]} mismatches={r[4]}\n")
    else:
        sys.stderr.write(
            f"[make_primer_bed] WARNING: primer not found in {args.ref} "
            f"(<= {args.max_mismatch} mismatches); empty BED written, nothing will be trimmed\n")


if __name__ == "__main__":
    main()
