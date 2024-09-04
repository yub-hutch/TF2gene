
## Install

`devtools::install_github('yub-hutch/TF2gene', auth_token = "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")`


## Utilities

### ATAC-seq peaks

- Summarize peaks in terms of distance to TSS & GC content: `summarize_consensus_peaks`
- Read bed file of real consensus peaks: `format_consensus_peaks`
- Get coordinaltes of null peaks: `get_control_peaks`
- Write control peaks to bed file: `write_control_peak_to_bed`

### Link ATAC-seq peaks to genes

- Calculate distance between peaks & genes: `calc_peak2gene_distance`


## Build a database of Cluster-Buster score on matched control peaks

### 1. Generate control peaks

```
null_peak = get_control_peaks(
  meta_gene = grch38,
  chr_lens = hg38_chr_lens,
  peak_length = 500,
  radius = 5e+05
)

write_peak_to_bed(null_peak, fbed = 'null_peaks.bed')

# bedtools getfasta -fi /fh/fast/sun_w/kenny_zhang/hg38_ref_genome.fa -bed null_peaks.bed -fo null_peaks.fa
```

### 2. Summarize control peaks

```
meta_null_peak = summarize_consensus_peaks(fasta = 'null_peaks.fa', ncores = 36)

saveRDS(meta_null_peak, file = 'meta_null_peak.rds')
```
