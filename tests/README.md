# Test data

```bash
nextflow run main.nf -profile test,docker
```

runs the whole pipeline on one small sample (`conf/test.config`: `--mode amplicon --rhinovirus`) and writes to `test_output/`. It takes a few minutes on a laptop.

## Data

`data/SRR40432756_random20k.fastq.gz` is a random subset of 20,000 reads of the public SRA run [SRR40432756](https://www.ncbi.nlm.nih.gov/sra/SRR40432756): a rhinovirus A51 amplicon library, single-end, with no host reads. It was drawn with the same seqtk the pipeline uses:

```bash
seqtk sample -s 42 SRR40432756.fastq.gz 20000 | gzip -n > SRR40432756_random20k.fastq.gz
```

The subset must be random: the reads in SRR40432756 are ordered by genome position, so its first 20,000 reads all fall in the 5'UTR and assemble nothing.

The run is distributed in SRA Lite format, where every base quality is replaced by Q30. The test therefore checks assembly, primer trimming and VP1 typing, but not quality filtering.

## Expected result

| Column | Value |
|---|---|
| `ref_tag` | `RV-A51` (OP342725.1) |
| `coverage_c1` | `100` |
| `mean_depth_c1` | `271.63` |
| `consensus_length` | `6705` |
| `pct_ns` / `num_ambiguous` | `0.0000` / `0` |
| `vp1_genotype` | `RV-A51`, 97.07% identity over 888 nt |

`expected/TEST_RVA51_OP342725.1_RV-A51_consensus_final.fa` is the consensus this run produces. With the same containers it is reproduced exactly:

```bash
cmp test_output/final_consensus_fasta/TEST_RVA51_OP342725.1_RV-A51_consensus_final.fa \
    tests/expected/TEST_RVA51_OP342725.1_RV-A51_consensus_final.fa
```

With a single sample, the consensus-vs-consensus check has nothing to compare, so `test_contamination.tsv` holds only its header.
