#!/usr/bin/env python3
"""
Build a NUCLEOTIDE VP1 typing database from a curated PICORNA genome FASTA.

The VP1 region of each rhinovirus genome (tags starting with 'RV-') is located
with a translated search (blastx) against an amino-acid VP1 reference set, then
the nucleotide region is extracted and extended to the full VP1 CDS span. Each
entry is labelled with the genome's CURATED tag (not the blastx hit), so the
labels are your curated truth; blastx is used only to locate the region.

Re-run this whenever the genome database changes -- e.g. after adding a newly
described genotype. Otherwise the two databases drift apart and the new
genotype is assembled but typed as its nearest neighbour:

    bin/build_vp1_nt_db.py \\
        --picorna assets/PICORNA_DB_v20280814_tag.fasta \\
        --vp1_aa  db/VP1_AA_164_annotated.fasta \\
        --out     db/VP1_NT_from_picorna.fasta

See "Adding a new genotype" in the README for the full checklist.

Requires blastx and makeblastdb (BLAST+) on PATH.
"""

import argparse
import os
import subprocess
import sys
import tempfile


def read_fasta(path):
    seqs, order = {}, []
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                header = line[1:].rstrip("\n")
                name = header.split()[0]
                buf = []
                order.append((name, header))
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs, order


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--picorna", required=True, help="curated PICORNA genome FASTA (headers: ACC<space>TAG)")
    ap.add_argument("--vp1_aa", required=True, help="amino-acid VP1 reference FASTA (used only to locate VP1)")
    ap.add_argument("--out", required=True, help="output nucleotide VP1 FASTA")
    ap.add_argument("--species_prefix", default="RV-", help="only extract from genomes whose tag starts with this (default: RV-)")
    ap.add_argument("--evalue", default="1e-10")
    ap.add_argument("--makeblastdb", action="store_true", help="also run makeblastdb -dbtype nucl on the output")
    args = ap.parse_args()

    for tool in ("blastx",) + (("makeblastdb",) if args.makeblastdb else ()):
        if subprocess.run(["which", tool], stdout=subprocess.DEVNULL).returncode != 0:
            sys.exit(f"ERROR: '{tool}' not found on PATH (install BLAST+).")

    genomes, order = read_fasta(args.picorna)
    # tag = second whitespace field of each header
    tags = {}
    with open(args.picorna) as fh:
        for line in fh:
            if line.startswith(">"):
                p = line[1:].split()
                tags[p[0]] = p[1] if len(p) > 1 else "NA"

    aa, _ = read_fasta(args.vp1_aa)
    slen = {k: len(v) for k, v in aa.items()}

    rv = [n for n in genomes if tags.get(n, "").startswith(args.species_prefix)]
    if not rv:
        sys.exit(f"ERROR: no genomes with tag prefix '{args.species_prefix}' found.")
    sys.stderr.write(f"[build_vp1_nt_db] {len(rv)} genomes with prefix '{args.species_prefix}'\n")

    with tempfile.TemporaryDirectory() as tmp:
        qfa = os.path.join(tmp, "rv.fa")
        with open(qfa, "w") as o:
            for n in rv:
                o.write(f">{n}\n{genomes[n]}\n")
        # blastx: locate VP1 (best HSP per genome)
        subprocess.run(["makeblastdb", "-in", args.vp1_aa, "-dbtype", "prot"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        hits = os.path.join(tmp, "hits.tsv")
        with open(hits, "w") as hf:
            subprocess.run(
                ["blastx", "-query", qfa, "-db", args.vp1_aa, "-max_target_seqs", "1",
                 "-max_hsps", "1", "-evalue", args.evalue,
                 "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send bitscore"],
                stdout=hf, stderr=subprocess.DEVNULL)

        best = {}
        for line in open(hits):
            f = line.rstrip("\n").split("\t")
            q, bs = f[0], float(f[8])
            if q not in best or bs > best[q][8]:
                best[q] = [f[0], f[1], float(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6]), int(f[7]), bs]

    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    n_written, no_hit = 0, []
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    with open(args.out, "w") as out:
        for n in rv:
            if n not in best:
                no_hit.append(n)
                continue
            _, ssid, pid, alen, qs, qe, ss, se, bs = best[n]
            g = genomes[n]
            L = len(g)
            sl = slen.get(ssid, se)
            # all picornavirus genomes are +sense; extend HSP to full VP1 CDS span
            if qs <= qe:
                lo = max(1, qs - (ss - 1) * 3)
                hi = min(L, qe + (sl - se) * 3)
                sub = g[lo - 1:hi]
            else:
                lo = max(1, qe - (sl - se) * 3)
                hi = min(L, qs + (ss - 1) * 3)
                sub = g[lo - 1:hi].translate(comp)[::-1]
            out.write(f">{n}_{tags[n]}\n{sub}\n")
            n_written += 1

    sys.stderr.write(f"[build_vp1_nt_db] wrote {n_written} VP1 nt sequences to {args.out}\n")
    if no_hit:
        sys.stderr.write(f"[build_vp1_nt_db] WARNING: no VP1 located for {len(no_hit)} genomes: {', '.join(no_hit)}\n")

    if args.makeblastdb:
        subprocess.run(["makeblastdb", "-in", args.out, "-dbtype", "nucl"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        sys.stderr.write(f"[build_vp1_nt_db] makeblastdb (nucl) done on {args.out}\n")


if __name__ == "__main__":
    main()
