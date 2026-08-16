# picorna-strm

**Reference-based assembly and VP1 genotyping of picornavirus genomes from short-read metagenomic, capture and amplicon sequencing.**

`picorna-strm` is **based on [revica-strm](https://github.com/greninger-lab/revica-strm)** (Greninger Lab, UW Virology), adapted and re-tuned for picornaviruses — small RNA virus genomes with extraordinary genetic diversity. It keeps revica-strm's two-pass consensus strategy and adds a curated picornavirus reference database, sequencing-mode presets (amplicon / capture / metagenomic), amplicon primer trimming, PCR-duplicate removal, VP1 genotyping, variant calling, and a number of correctness fixes.

This is an independent fork: the pipeline is expected to keep diverging from upstream as it is further adapted to picornavirus biology.

---

## What it does

For each sample, the pipeline runs:

1. **Quality filtering** — `fastp` (adapter and quality trimming). Always runs.
2. *(optional)* **Host removal** — Kraken2 against a host database **you supply** (`--run_kraken2 --kraken2_db`); none is bundled. See `docs/making_kraken2_human_db.md`.
3. *(optional)* **Subsampling** — `seqtk`, 6M reads by default.
4. **Reference detection** — maps reads against the whole picornavirus database and keeps every reference passing coverage/depth thresholds. **More than one reference can be selected per sample**, so co-infections yield one assembly per organism.
5. **First consensus** — reads are mapped to each selected reference and `iVar` builds a rough consensus. Depending on `--mode`, the alignment is primer-trimmed (amplicon) or deduplicated (capture) first.
6. **Final consensus** — reads are re-mapped against that first consensus (now a reference of the same organism that was actually sequenced), the same mode-specific correction is applied again, and `iVar` builds the delivered assembly.
7. **Per-position VCF** — `bcftools mpileup` on that final alignment against the sample's own consensus, recording intra-host diversity at every covered position.
8. *(optional)* **VP1 genotyping** — each final consensus is typed by `blastn` against a nucleotide VP1 database.

VP1 typing is deliberately **decoupled** from reference selection: assembly still uses whole-genome mapping, and VP1 provides the authoritative genotype call on top of the finished consensus. Rhinovirus types are defined by VP1 divergence, so VP1 can resolve genotype pairs that whole-genome mapping cannot separate.

## Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/install.html)
- [Docker](https://docs.docker.com/desktop/) (all bioinformatics tools ship in containers)
- Python 3.8+ (for the wrapper script)

> **Do not run the pipeline from a cloud-synced folder** (OneDrive, iCloud Drive, Dropbox). Docker cannot read macOS File Provider mounts: input files appear as 0 bytes inside the container and tools fail with confusing errors such as `igzip: Error invalid gzip header`. Keep the pipeline, the FASTQ files and the work directory on local disk.

## Quick start

```bash
git clone <your-repo-url>
cd picorna-strm
chmod +x picorna-strm bin/*.py     # only needed if the exec bit was lost

./picorna-strm /path/to/fastq_dir -profile docker --output results
```

The first non-flag argument is the FASTQ directory; it must come **before** any option that takes a value.

Typical runs:

```bash
# amplicon: forward primer trimmed from consensus and variants
./picorna-strm /path/to/amplicon_fastqs -profile docker \
    --mode amplicon --output results_amplicon --rhinovirus

# hybridization capture: PCR duplicates removed
./picorna-strm /path/to/capture_fastqs -profile docker \
    --mode capture --output results_capture --rhinovirus --sample false

# shotgun metagenomics: no primers, few PCR cycles -> neither
./picorna-strm /path/to/metagenomic_fastqs -profile docker \
    --mode metagenomic --output results_meta --rhinovirus --sample false
```

Run `./picorna-strm -h` for all options.

### FASTQ naming

Sample names must be unique **before the first underscore**:

- correct: `sample1_R1.fastq.gz`, `sample1_R2.fastq.gz`
- wrong: `sample_1_R1.fastq.gz`, `sample_2_R1.fastq.gz`

If filenames contain no underscore, the same logic applies to the first `.`.

## Sequencing modes

`--mode` tells the pipeline how the library was made; everything library-specific follows from it.

| `--mode` | Primer trimming | Duplicate removal | Use for |
|---|---|---|---|
| `amplicon` | yes | **no** | targeted amplicon sequencing |
| `capture` | no | yes | hybridization capture / enrichment |
| `metagenomic` | no | no | shotgun metagenomics (default) |

The pairing is not arbitrary. Amplicon reads share start coordinates by design, so deduplicating them would collapse genuine molecules and destroy the depth; conversely their primer-derived bases reflect the primer rather than the sample and must be removed. Capture libraries have no primers to trim but do accumulate PCR duplicates that inflate apparent depth. Shotgun metagenomics has neither problem, since there are no primers and few PCR cycles.

Both operations are applied to **each of the two alignment passes**, so the first consensus, the final consensus and the variant calls are all derived from the same corrected alignment.

### Primer trimming

Trimming uses `ivar trim` (iVar 1.4.4; 1.4.1-1.4.3 fixed several `trim` bugs, so older versions are not recommended). Trimming is by **coordinate**, not by matching the primer sequence in the read, so reads carrying only part of the primer are handled correctly. The coordinates are located per task with `bin/make_primer_bed.py`, which searches the reference of that particular alignment for `--primer_fwd` (IUPAC codes supported, up to `--primer_max_mismatch` mismatches, both strands). Doing it per task matters because the primer sits at a different coordinate in every genotype, and at a different coordinate again in the sample's own first consensus.

The default primer is the conserved 5'UTR forward primer `TCCTCCGGYCCCTGAATGYGGCTAA`. The reverse primer anneals to the genomic poly(A) tail, so there is nothing to trim at the 3' end. If the primer is not found in a given reference, an empty BED is written and the alignment passes through untrimmed rather than failing.

Reads that contain no primer are kept (`ivar trim -e`), so coverage away from the primer-binding site is unaffected.

## Per-position VCF (intra-host diversity)

With `--call_variants` (on by default), `bcftools mpileup` runs on the **same alignment the final consensus was built from** — primer-trimmed and/or deduplicated to match — using that sample's own first consensus as the reference. One VCF per consensus genome is written to `<output>/vcf/`.

This is deliberately *not* a variant caller: `bcftools call` is not run, so **every covered position** is reported with its allelic depths (`FORMAT/AD`) and total depth (`FORMAT/DP`). The file is therefore a complete record of **intra-host, intra-genotype diversity**, from which allele frequencies can be recomputed at any site downstream — including sites that no caller would have flagged. It reproduces:

```bash
bcftools mpileup -f <consensus> -a FORMAT/AD,FORMAT/DP \
    --max-depth 1000000 --min-BQ 20 --min-MQ 20 -Ov -o <sample>.vcf <bam>
```

Tunable with `--vcf_min_bq`, `--vcf_min_mq` and `--vcf_max_depth`.

Because the reference is the sample's own consensus, these are differences **within** the sample, not differences from a database reference. In `--mode amplicon` the primer-binding site is absent from the alignment, so no primer-derived position ever appears in the VCF.

### Optional: called-variant table

`--ivar_variants` additionally runs `ivar variants` on the same alignment, producing `*.variants.tsv` in `final_files/variants/` — a filtered table of variants with per-allele depth and frequency. Thresholds: `--ivar_var_q`, `--ivar_var_t`, `--ivar_var_m`. The VCF itself comes from `bcftools mpileup` above; iVar is used here only for its allele-frequency table.

## Output

Inside `--output`:

| Path | Contents |
|---|---|
| `final_files/final_assemblies/` | final consensus genomes (one per sample per detected organism) |
| `final_files/align_to_db/` | coverage of every database reference — useful to audit what was and was not detected |
| `final_files/align_to_selected_ref/` | BAMs against the selected reference(s) |
| `final_files/align_to_consensus/` | first consensus and BAMs of the second pass |
| `vcf/` | per-position VCF per consensus genome (intra-host diversity) |
| `final_files/variants/` | called-variant table, with `--ivar_variants` |
| `final_files/vp1_genotyping/` | per-consensus VP1 typing (with `--rhinovirus`) |
| `primer_trim/` | primer BED located for each alignment (`--mode amplicon`) |
| `final_files/SRA/` | mapped reads as FASTQ, ready for submission |
| `<run_name>_summary.tsv` | reads, coverage, depth, consensus length, %N per assembly |
| `<run_name>_vp1_genotypes.tsv` | consolidated VP1 genotype calls |
| `<run_name>_multiqc.html` | MultiQC report |
| `fail/` | references and samples that did not pass thresholds |

## Databases

### Genome database (assembly)

`assets/PICORNA_DB_v20280814_tag.fasta` — curated picornavirus genomes (enterovirus A–D, rhinovirus A/B/C, coxsackievirus, echovirus, poliovirus). Select a different one with `--refs`.

FASTA headers must follow:

```
>ACCESSION<space>TAG[<space>DESCRIPTION]
```

`TAG` must be unique per genotype/segment — the pipeline keeps **one reference per tag**, so a tag shared across genotypes collapses them into a single call. The description field is optional.

```
>OP342726.1 RV-A49
>K02121.1 RV-B14
```

### VP1 database (genotyping)

`db/VP1_NT_from_picorna.fasta` — nucleotide VP1 sequences extracted from the genome database above, each labelled with its curated tag. Typing is done with `blastn`; the BLAST index is built inside the task, so no index files are stored in the repository.

**Regenerate this whenever the genome database changes**, or the two will drift apart:

```bash
bin/build_vp1_nt_db.py \
    --picorna assets/PICORNA_DB_v20280814_tag.fasta \
    --vp1_aa  db/VP1_AA_164_annotated.fasta \
    --out     db/VP1_NT_from_picorna.fasta
```

`db/VP1_AA_164_annotated.fasta` (amino-acid VP1 references) is used **only** by that builder, to locate the VP1 region within each genome. Genotype labels always come from the curated genome tags, never from the amino-acid hit. The pipeline itself never reads it.

---

## Adding a new genotype

New rhinovirus (and enterovirus) types keep being described. The genome database and the VP1 typing database **must be updated together** — the VP1 database is *derived* from the genome database, so editing one without rebuilding the other silently leaves them out of sync, and the new genotype will be assembled but typed as its nearest neighbour.

### 1. Add the genome to the reference database

Append the new genome to `assets/PICORNA_DB_v20280814_tag.fasta` using the required header format:

```
>ACCESSION<space>TAG[<space>DESCRIPTION]
```

```
>PV717813.1 RV-A44
```

The `TAG` must be **unique to the genotype** — the pipeline keeps only one reference per tag, so reusing a tag for two genotypes collapses them into a single call.

### 2. Check database integrity

Two failure modes are silent and worth checking every time:

```bash
DB=assets/PICORNA_DB_v20280814_tag.fasta

# (a) duplicate accessions -- samtools faidx SILENTLY DROPS one of them
grep "^>" $DB | awk '{print $1}' | sort | uniq -d

# (b) malformed headers (must have at least accession + tag)
grep "^>" $DB | awk 'NF<2'

# (c) confirm every sequence is indexable
samtools faidx $DB && wc -l < ${DB}.fai   # must equal the number of ">" lines
```

Both commands must print nothing. A duplicate accession is the most damaging error: the pipeline will extract the wrong sequence for one of the two entries without any warning.

### 3. Rebuild the VP1 typing database

```bash
bin/build_vp1_nt_db.py \
    --picorna assets/PICORNA_DB_v20280814_tag.fasta \
    --vp1_aa  db/VP1_AA_164_annotated.fasta \
    --out     db/VP1_NT_from_picorna.fasta
```

The builder locates VP1 in every genome tagged `RV-*` (translated search against the amino-acid references), extracts the nucleotide region, extends it to the full CDS span and labels it with the curated tag. It reports how many sequences were written and warns about any genome where VP1 could not be located.

> A genotype whose VP1 is too divergent to be found by the amino-acid references will be skipped with a warning. If that happens, add a closer amino-acid VP1 reference to `db/VP1_AA_164_annotated.fasta` and rebuild. To extract VP1 from a non-rhinovirus species, pass e.g. `--species_prefix EV-`.

### 4. Verify the two databases are in sync

```bash
grep "^>" assets/PICORNA_DB_v20280814_tag.fasta \
  | awk '$2 ~ /^RV-/ {print substr($1,2)"_"$2}' | sort > /tmp/genomes.txt
grep "^>" db/VP1_NT_from_picorna.fasta | sed 's/^>//' | sort > /tmp/vp1.txt
diff /tmp/genomes.txt /tmp/vp1.txt && echo "in sync"
```

`diff` must report no differences: every rhinovirus genome should have exactly one VP1 entry.

### 5. Sanity-check the new genotype

Confirm the new entry types as itself, and see what it is closest to:

```bash
makeblastdb -in db/VP1_NT_from_picorna.fasta -dbtype nucl -out /tmp/vp1db
samtools faidx assets/PICORNA_DB_v20280814_tag.fasta <ACCESSION> > /tmp/new.fa
blastn -query /tmp/new.fa -db /tmp/vp1db -max_target_seqs 5 -max_hsps 1 \
       -perc_identity 72 -outfmt "6 sseqid pident length bitscore"
```

The top hit must be the new genotype itself at ~100%. If the **second** hit is a *different* genotype above ~97% identity, the two are hard to separate by VP1 and calls between them will be flagged `ambiguous` at run time — worth noting, and worth double-checking that the new entry is not simply a mislabelled duplicate of an existing one.

Finally, commit both files together so the databases never drift apart in version control:

```bash
git add assets/PICORNA_DB_v20280814_tag.fasta db/VP1_NT_from_picorna.fasta
git commit -m "Add RV-Axx (<accession>); rebuild VP1 typing database"
```

## Key parameters

| Parameter | Default | Notes |
|---|---|---|
| `--refs` | bundled PICORNA db | genome database for reference detection |
| `--sample` | `6000000` | subsampling; use `false` to disable (recommended for unbalanced co-infections) |
| `--ref_min_cov` / `--ref_min_depth` | `30` / `3` | reference selection thresholds. For high-coverage capture data, `--ref_min_cov 60` suppresses spurious relatives of genotypes missing from the database |
| `--ivar_fin_m` | `5` | minimum depth for a final consensus base. At ~10x mean depth, a value of 10 masks ~47% of a minor co-infecting genome as N versus ~8% at 5 |
| `--mode` | `metagenomic` | `amplicon` \| `capture` \| `metagenomic` — see above |
| `--primer_fwd` | 5'UTR primer | forward primer for `--mode amplicon` (IUPAC allowed) |
| `--call_variants` | on | per-position VCF (`vcf/`) against the final consensus |
| `--rhinovirus` | off | enable VP1 genotyping |
| `--run_kraken2` + `--kraken2_db` | off | remove host reads before assembly (see `docs/making_kraken2_human_db.md`) |
| `--vp1_min_identity` / `--vp1_min_aln_len` | `72` / `300` | minimum blastn identity (%) and alignment length (nt) for a VP1 hit |
| `--fastp_container` | biocontainers 0.23.2 | override if the default image misbehaves on your platform |

### Interpreting VP1 calls

The winner is the hit with the highest **bitscore** (which combines identity and alignment length), not simply the highest identity. The 72% identity floor is a noise filter, **not a taxonomic threshold**: a genotype absent from the database will be reported as its nearest relative rather than as a novel type. Treat identities below ~90% as candidate unrepresented types. When the runner-up genotype falls within `--vp1_ambiguous_margin` (1% by default), the call is flagged `ambiguous:<genotype>(<identity>%)` in the `vp1_note` column.

A `no_hit` result is not a failure — non-rhinovirus consensuses (enterovirus, echovirus, …) simply have no VP1 entry in this database. VP1 typing never gates assembly.

## Known limitations

- Two strains of the **same** genotype in one sample collapse into a single assembly (one reference per tag).
- Only the **forward** primer is trimmed in `--mode amplicon`, which matches a design whose reverse primer anneals to the poly(A) tail. A multi-amplicon (tiling) scheme would need a full primer BED instead.
- A genotype present in the sample but **absent from the genome database** scatters its reads across related genotypes, several of which may pass the selection thresholds with partial (~35–55%) coverage. Genuine detections typically show ~99–100% coverage.

## Credits

This pipeline is a modified fork of [**revica-strm**](https://github.com/greninger-lab/revica-strm) by Eli Piliper, Jaydee Sereewit, Stephanie Goya and Alex L. Greninger (UW Medicine, Department of Laboratory Medicine and Pathology, University of Washington), which is itself derived from [REVICA](https://github.com/greninger-lab/revica).

Licensed under the GNU General Public License, following upstream.
