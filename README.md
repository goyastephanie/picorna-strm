# picorna-strm

**Reference-based assembly and VP1 genotyping of picornavirus genomes from short-read metagenomic, capture and amplicon sequencing.**

`picorna-strm` is **based on [revica-strm](https://github.com/greninger-lab/revica-strm)** (Greninger Lab, UW Virology), adapted and re-tuned for picornaviruses — small RNA virus genomes with extraordinary genetic diversity. It keeps revica-strm's two-pass consensus strategy and adds a curated picornavirus reference database, sequencing-mode presets (amplicon / capture / metagenomic), amplicon primer trimming, VP1 genotyping, per-position variant calling, a consensus-vs-consensus contamination and duplicate-assembly check, a software version registry, and a number of correctness fixes.

This is an independent fork: the pipeline is expected to keep diverging from upstream as it is further adapted to picornavirus biology.

---

## What it does

For each sample, the pipeline runs:

1. **Quality filtering** — `fastp` (adapter and quality trimming). Always runs.
2. *(optional)* **Host removal** — Kraken2 against a host database **you supply** (`--run_kraken2 --kraken2_db`); none is bundled. See `docs/making_kraken2_human_db.md`.
3. *(optional)* **Subsampling** — `seqtk`, 6M reads by default.
4. **Reference detection** — maps reads against the whole picornavirus database and keeps every reference passing coverage/depth thresholds. **More than one reference can be selected per sample**, so co-infections yield one assembly per organism.
5. **First consensus** — reads are mapped to each selected reference and `iVar` builds a rough consensus. In `--mode amplicon` the alignment is primer-trimmed first.
6. **Final consensus** — reads are re-mapped against that first consensus (now a reference of the same organism that was actually sequenced), the same mode-specific correction is applied again, and `iVar` builds the delivered assembly.
7. **Per-position VCF** — `bcftools mpileup` on that final alignment against the sample's own consensus, recording intra-host diversity at every covered position.
8. *(optional)* **VP1 genotyping** — each final consensus is typed by `blastn` against a nucleotide VP1 database.
9. **Consensus-vs-consensus check** — all final consensuses of the run are compared against each other: near-identical pairs from *different* samples are flagged as possible contamination, and near-identical pairs from the *same* sample as one virus assembled against two genotype references.

VP1 typing is deliberately **decoupled** from reference selection: assembly still uses whole-genome mapping, and VP1 provides the authoritative genotype call on top of the finished consensus. Rhinovirus types are defined by VP1 divergence, so VP1 can resolve genotype pairs that whole-genome mapping cannot separate.

## Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/install.html)
- [Docker](https://docs.docker.com/desktop/) (all bioinformatics tools ship in containers)
- Python 3.8+ (for the wrapper script)

> **Do not run the pipeline from a cloud-synced folder** (OneDrive, iCloud Drive, Dropbox). Docker cannot read macOS File Provider mounts: input files appear as 0 bytes inside the container and tools fail with confusing errors such as `igzip: Error invalid gzip header`. Keep the pipeline, the FASTQ files and the work directory on local disk.

## Quick start

```bash
git clone https://github.com/goyastephanie/picorna-strm.git
cd picorna-strm
chmod +x picorna-strm bin/*.py

./picorna-strm /path/to/fastq_dir -profile docker --output results
```

The first non-flag argument is the FASTQ directory; it must come **before** any option that takes a value.

Typical runs:

```bash
# amplicon sequencing of rhinovirus genomes: forward primer trimmed from consensus and variants
./picorna-strm /path/to/amplicon_fastqs -profile docker \
    --mode amplicon --output results_amplicon --rhinovirus

# hybridization capture of enterovirus (currently processed like metagenomic)
./picorna-strm /path/to/capture_fastqs -profile docker \
    --mode capture --output results_capture

# shotgun metagenomics of rhinovirus: no primers, few PCR cycles -> neither
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

| `--mode` | Primer trimming | Use for |
|---|---|---|
| `amplicon` | yes | targeted amplicon sequencing |
| `capture` | no | hybridization capture / enrichment |
| `metagenomic` | no | shotgun metagenomics (default) |

Amplicon reads carry primer-derived bases that reflect the primer pool rather than the sample, so they must be removed; deduplication, on the other hand, would be actively harmful here, because amplicon reads share start coordinates by design.

`capture` currently behaves exactly like `metagenomic`. It exists as its own mode so that capture-specific handling can be added later without changing anyone's command line. Coordinate-based deduplication was evaluated and removed: on a 7 kb genome sequenced single-end, read start positions saturate, and duplicate marking discarded ~97% of real reads in test data.

Primer trimming is applied to **each of the two alignment passes**, so the first consensus, the final consensus and the variant calls all come from the same corrected alignment.

### Primer trimming

Trimming uses `ivar trim` (iVar 1.4.4, pinned in the container). Trimming is by **coordinate**, not by matching the primer sequence in the read, so reads carrying only part of the primer are handled correctly. The coordinates are located per task with `bin/make_primer_bed.py`, which searches the reference of that particular alignment for `--primer_fwd` (IUPAC codes supported, up to `--primer_max_mismatch` mismatches, both strands). The search is repeated for each alignment pass because the two references are different sequences.

Against the **database reference** the primer sits at its genomic coordinate — around position 300-450 of the 5'UTR, and different in every genotype — so it is found and trimmed. This is the primary trim, and doing it here matters: the first consensus becomes the reference for the second alignment and for the VCF, so primer-derived bases must not reach it. Left untrimmed they are clearly visible, including IUPAC ambiguity codes at the primer's degenerate positions (`TCCTCCGG**Y**CCCTGAATG**Y**GGCTAA`) — the sample has one base there, the primer pool has both.

Against the **sample's own first consensus** the search acts as a safety net. If pass 1 worked, the consensus begins just downstream of the primer, the BED comes back empty and nothing is trimmed — correct, because reads still carrying the primer have those bases hanging off the 5' end and the aligner soft-clips them anyway. If any primer did survive into the first consensus it is caught here; since quality trimming leaves a variable amount of primer on each read, that residual is often *truncated*, so the search also matches a primer overhanging the end of the reference (the `--min_overlap` option of `bin/make_primer_bed.py`, 10 bases; not exposed as a pipeline flag). Both paths give a primer-free final consensus and VCF.

The default primer is the conserved 5'UTR forward primer `TCCTCCGGYCCCTGAATGYGGCTAA`. The reverse primer anneals to the genomic poly(A) tail, so there is nothing to trim at the 3' end. If the primer is not found in a given reference, an empty BED is written and the alignment passes through untrimmed rather than failing.

Reads that contain no primer are kept (`ivar trim -e`), so coverage away from the primer-binding site is unaffected.

## Per-position VCF (intra-host diversity)

With `--call_variants` (on by default), `bcftools mpileup` runs on the **same alignment the final consensus was built from** — primer-trimmed to match, in `--mode amplicon` — using that sample's own first consensus as the reference. One VCF per consensus genome is written to `<output>/vcf/`.

This is deliberately *not* a variant caller: `bcftools call` is not run, so **every covered position** is reported with its allelic depths (`FORMAT/AD`) and total depth (`FORMAT/DP`). The file is therefore a complete record of **intra-host, intra-genotype diversity**, from which allele frequencies can be recomputed at any site downstream — including sites that no caller would have flagged. It reproduces:

```bash
bcftools mpileup -f <consensus> -a FORMAT/AD,FORMAT/DP \
    --max-depth 1000000 --min-BQ 20 --min-MQ 20 -Ov -o <sample>.vcf <bam>
```

Tunable with `--vcf_min_bq`, `--vcf_min_mq` and `--vcf_max_depth`.

Because the reference is the sample's own consensus, these are differences **within** the sample, not differences from a database reference. In `--mode amplicon` the primer-binding site is absent from the alignment, so no primer-derived position ever appears in the VCF.

### Optional: called-variant table

`--ivar_variants` additionally runs `ivar variants` on the same alignment, producing `*.variants.tsv` in `variants/` — a filtered table of variants with per-allele depth and frequency. Thresholds: `--ivar_var_q`, `--ivar_var_t`, `--ivar_var_m`. The VCF itself comes from `bcftools mpileup` above; iVar is used here only for its allele-frequency table.

## Output

Inside `--output`:

| Path | Contents |
|---|---|
| `final_consensus_fasta/` | **the delivered consensus genomes**, one per sample per detected organism |
| `final_consensus_assemblies/` | the alignment each of those consensuses was called from: the second-pass BAM, its `.bai`, and the **first consensus** in FASTA, which is that BAM's reference and the reference of every VCF. In `--mode amplicon` the BAM published here is the primer-trimmed one, because that is the alignment the consensus and the VCF actually come from; the untrimmed version is not published |
| `vcf/` | per-position VCF per consensus genome (intra-host diversity), in the coordinates of the first consensus |
| `align_to_db/` | `*_covstats.tsv`: coverage of every database reference — the evidence behind reference selection. The BAM only with `--save_db_bam` |
| `align_to_selected_ref/` | first-pass BAM, `.bai` and the selected database reference in FASTA. Only with `--save_ref_bam` |
| `variants/` | called-variant table, with `--ivar_variants` |
| `primer_trim/` | primer BED located for each alignment (`--mode amplicon`) |
| `SRA/` | mapped reads as FASTQ, only with `--save_sra_fastq` |
| `fastp/json/` | fastp JSON per sample: the machine-readable QC record (Q20/Q30, duplication, adapter content, per-cycle quality) and the only way to regenerate the MultiQC report after `work/` is deleted |
| `kraken2/` | Kraken2 reports, with `--run_kraken2` |
| `seqtk_sample/` | subsampled FASTQs, with `--save_sample_reads` |
| `<run_name>_summary.tsv` | 19 columns per assembly: raw/trimmed reads, mapped reads, coverage and mean depth for both alignment passes, consensus length, Ns and ambiguity codes |
| `<run_name>_vp1_genotypes.tsv` | consolidated VP1 genotype calls |
| `<run_name>_contamination.tsv` | all-vs-all comparison of the run's consensuses; flags possible contamination and duplicate assemblies |
| `<run_name>_multiqc.html` | MultiQC report |
| `fail/` | references and samples that did not pass thresholds |
| `pipeline_info/software_versions.yml` | container image per process, Nextflow build, pipeline revision, database MD5s |
| `params.json` | every parameter the run actually used |

The two FASTA files per assembly are **not** duplicates. `final_consensus_fasta/*_consensus_final.fa` is the delivered genome, called from the second alignment pass. `final_consensus_assemblies/*_consensus1.fa` is the intermediate scaffold: it is the reference the second-pass BAM is aligned to and the reference the VCF coordinates and REF alleles refer to, so deleting it makes the BAMs and VCFs uninterpretable. That is why it lives next to them rather than next to the deliverable.

### Disk footprint

Everything is published as a **copy**, and every optional output is gated with `saveAs` returning `null` rather than with `enabled`. 

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

**Regenerate this whenever the genomes database changes**, or the two will drift apart:

```bash
bin/build_vp1_nt_db.py \
    --picorna assets/PICORNA_DB_v20280814_tag.fasta \
    --vp1_aa  db/VP1_AA_164_annotated.fasta \
    --out     db/VP1_NT_from_picorna.fasta
```

`db/VP1_AA_164_annotated.fasta` (amino-acid VP1 references) is used **only** by that builder, to locate the VP1 region within each genome. Genotype labels always come from the curated genome tags, never from the amino-acid hit. The pipeline itself never reads it.

---

## Adding a new genotype

New rhinovirus and enterovirus types keep being described. The genome database and the rhinovirus VP1 typing database **must be updated together** — the VP1 database is *derived* from the genome database, so editing one without rebuilding the other silently leaves them out of sync, and the new genotype will be assembled but typed as its nearest neighbour.

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
| `--sample` | `6000000` | subsampling after fastp, fixed seed; use `false` to disable (recommended for unbalanced co-infections). Capped at `--amplicon_max_sample` in `--mode amplicon` |
| `--ref_min_cov` / `--ref_min_depth` | `30` / `3` | reference selection thresholds. For high-coverage capture data, `--ref_min_cov 60` suppresses spurious relatives of genotypes missing from the database |
| `--ivar_fin_m` | `5` | minimum depth for a final consensus base. At ~10x mean depth, a value of 10 masks ~47% of a minor co-infecting genome as N versus ~8% at 5 |
| `--ivar_init_t` / `--ivar_fin_t` | `0.4` / `0.6` | frequency a base must reach to be called on its own. A **higher** value produces **more** IUPAC ambiguity codes. Low for the scaffold consensus (it only needs to remove the divergence from a distant reference), higher for the final one (a 55/45 site is written as an ambiguity code, recording intra-sample diversity) |
| `--mpileup_max_depth` / `--mpileup_max_depth_final` | `100` / `500` | pileup depth cap for the first and the final consensus; `0` for unlimited |
| `--min_trimmed_reads` | `1000` | samples with fewer reads after `fastp` are dropped (a `WARN` line names each one) |
| `--mode` | `metagenomic` | `amplicon` \| `capture` \| `metagenomic` — see above |
| `--primer_fwd` | 5'UTR primer | forward primer for `--mode amplicon` (IUPAC allowed) |
| `--primer_max_mismatch` | `2` | mismatches tolerated when locating the primer in a reference |
| `--check_contamination` | on | all-vs-all comparison of the run's consensuses |
| `--contam_flag_identity` | `99.9` | flag two **different** samples above this identity as possible contamination |
| `--contam_same_genome_identity` | `99.0` | flag two assemblies of **one** sample above this identity as the same virus |
| `--contam_min_comparable` | `500` | minimum unambiguous overlap needed to judge a pair |
| `--contam_intraspecies_only` | on | compare only consensuses whose tag prefix (RV-A, RV-C, EV-B, …) matches |
| `--call_variants` | on | per-position VCF (`vcf/`), in the coordinates of the **first** consensus, which is the reference of the alignment it is built from |
| `--rhinovirus` | off | enable VP1 genotyping |
| `--run_kraken2` + `--kraken2_db` | off | remove host reads before assembly (see `docs/making_kraken2_human_db.md`) |
| `--vp1_min_identity` / `--vp1_min_aln_len` | `72` / `300` | **reporting** floor: minimum identity (%) and unambiguous comparable positions for a VP1 hit to be shown at all |
| `--vp1_type_identity` / `--vp1_type_identity_b` | `87.0` / `88.0` | **type assignment** threshold (RV-A and RV-C / RV-B). Below it the call is `unassigned` with the nearest relative named |
| `--amplicon_max_sample` | `3000000` | in `--mode amplicon`, cap on `--sample` |
| `--fastp_container` | biocontainers 0.23.2 | override if the default image misbehaves on your platform |

### Interpreting VP1 calls

Two thresholds do two different jobs. `--vp1_min_identity` (72%) is a **reporting** floor: below it the hit is noise and nothing is shown. `--vp1_type_identity` (87%, and `--vp1_type_identity_b` 88% for RV-B) is the **type assignment** threshold: [McIntyre et al. 2013](https://www.microbiologyresearch.org/content/journal/jgv/10.1099/vir.0.053686-0) place the inter-type/intra-type boundary at 13% VP1 nucleotide divergence for RV-A and RV-C and 12% for RV-B. A consensus that clears the reporting floor but not the assignment threshold is reported as `vp1_genotype = unassigned`, with its nearest relative and identity still in `vp1_subject` and `vp1_identity` and named in `vp1_note` — the right answer both for a genuinely novel type and for one simply missing from the database.

Among the entries clearing the reporting floor, the winner is the one with the highest **summed bitscore** across its non-overlapping HSPs, not simply the highest identity. The 72% identity floor is a noise filter, **not a taxonomic threshold**: a genotype absent from the database is reported as its nearest relative rather than as a novel type. Treat identities below ~90% as candidate unrepresented types. When the runner-up genotype falls within `--vp1_ambiguous_margin` (1% by default), the call is flagged `ambiguous:<genotype>(<identity>%)` in the `vp1_note` column.

`vp1_identity` and `vp1_aln_len` do not come straight from blast. blastn is a **local** aligner: it returns HSPs, and a run of Ns long enough to trip its X-drop cuts the alignment in two, after which blast re-seeds on the other side. Measured on this database with a real 849 nt VP1 and N runs spaced every 200 bases, runs of up to 25 Ns are absorbed into a single full-length HSP, while runs of 30 or more break it into ~200 nt fragments. Reporting only the best fragment would then fail the 300 nt floor for a consensus that is only 10% N. So all non-overlapping HSPs of an entry are summed, and identity is computed **only over positions where both sequences call an unambiguous A/C/G/T** — the same treatment used by the cross-contamination check, and for the same reason. `vp1_aln_len` is therefore the number of comparable unambiguous positions, not blast's raw alignment length.

The search runs with `-task blastn` (word size 11) rather than the default megablast (word size 28). Megablast loses the hit above roughly 13% divergence, which puts the 72% identity floor out of its reach: a divergent or unrepresented genotype came back as `no_hit` instead of as its nearest relative. Verified on all 184 database VP1s used as queries: identical genotype calls under both, with no new ambiguous flags.

A `no_hit` result is not a failure — non-rhinovirus consensuses (enterovirus, echovirus, …) simply have no VP1 entry in this database. VP1 typing never gates assembly.

## Known limitations

- Two strains of the **same** genotype in one sample collapse into a single assembly (one reference per tag).
- Conversely, one virus can be assembled **twice** against two different but closely related genotype references of its species, looking like a co-infection. The consensus-vs-consensus check flags this as `possible_same_genome`; it is detected and reported, not prevented.
- Only the **forward** primer is trimmed in `--mode amplicon`, which matches a design whose reverse primer anneals to the poly(A) tail. A multi-amplicon (tiling) scheme would need a full primer BED instead.
- The consensus and the VCF for intra-sample study apply **different mapping-quality filters**: `iVar` builds the consensus from a pileup with no MAPQ floor (the iVar convention: `--min-BQ 0` in mpileup, quality filtering delegated to `ivar -q` including a minimum depth of 5X), while `bcftools mpileup` uses `--min-MQ 20` and a minimum depth 20X. On a 7 kb non-repetitive genome the two agree in practice, but at a repetitive locus the VCF can be blank where the consensus called a base.
- A genotype present in the sample but **absent from the genome database** scatters its reads across related genotypes, several of which may pass the selection thresholds with partial (~35–55%) coverage. Genuine detections typically show ~99–100% coverage.

## Credits

This pipeline is a modified fork of [**revica-strm**](https://github.com/greninger-lab/revica-strm) by Eli Piliper, Jaydee Sereewit, Stephanie Goya and Alex L. Greninger (UW Medicine, Department of Laboratory Medicine and Pathology, University of Washington), which is itself derived from [REVICA](https://github.com/greninger-lab/revica).

Licensed under the GNU General Public License, following upstream.
