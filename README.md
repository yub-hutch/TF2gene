
## 1. Install

`devtools::install_github('yub-hutch/TF2gene', auth_token = "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")`


## 2. Utilities

### ATAC-seq peaks

- Summarize peaks in terms of distance to TSS & GC content: `summarize_consensus_peaks`
- Read bed file of real consensus peaks: `format_consensus_peaks`
- Get coordinaltes of null peaks: `get_control_peaks`
- Write control peaks to bed file: `write_control_peak_to_bed`

### Link ATAC-seq peaks to genes

- Calculate distance between peaks & genes: `calc_peak2gene_distance`


## 3. Build a database of Cluster-Buster score on matched control peaks

### 1. Generate control peaks

```
null_peak = get_control_peaks(
  meta_gene = grch38,
  chr_lens = hg38_chr_lens,
  peak_length = 500,
  radius = 5e+05
)

write_control_peak_to_bed(null_peak, fbed = 'null_peaks.bed')
```

```
bedtools getfasta -fi /fh/fast/sun_w/kenny_zhang/hg38_ref_genome.fa -bed null_peaks.bed -fo null_peaks.fa
```

5,807,331 control peaks in total.

### 2. Summarize control peaks

```
meta_null_peak = summarize_consensus_peaks(fasta = 'null_peaks.fa', ncores = 36)
```

```
ggplot(meta_null_peak, aes(dist, gc_content)) +
  geom_hex(bins = 100) +
  labs(x = 'Distance to nearest TSS', y = 'GC content') +
  stat_cor(method = 'spearman', label.x.npc = 'center') +
  guides(fill = 'none') +
  theme_pubr()
```

![image](https://github.com/user-attachments/assets/fc180974-f8f8-4e0e-960a-002805a9a711)

### 3. Run Cluster-Buster on control peaks

Do motif by motif.

```
while IFS= read -r motif; do
  echo "$motif" > "feather/$motif.lst"
done < /fh/fast/sun_w/kenny_zhang/v10nr_clust_public/motifs.lst
```

```
#!/bin/bash

head -n 2 /fh/fast/sun_w/kenny_zhang/v10nr_clust_public/motifs.lst | while IFS= read -r motif; do
  sbatch --job-name=CB_contol \
         --nodes=1 \
         --ntasks-per-node=1 \
         --cpus-per-task=1 \
         --time=1:00:00 \
         --mem=10G \
         --output=/dev/null \
         --error=/dev/null \
         --wrap='python3 /fh/fast/sun_w/kenny_zhang/create_cisTarget_databases/create_cistarget_motif_databases.py \
                         -f null_peaks.fa \
                         -M /fh/fast/sun_w/kenny_zhang/v10nr_clust_public/singletons \
                         -m feather/$motif.lst \
                         -o feather/$motif
                         -t 1'
done
```

## 4. Calculate P-value for Cluster-buster score

To sample matched controls of ATAC-seq consensus peaks, generate weights for control peaks in 3 steps:
1. fit (distance to nearest TSS, GC content) joint density $D$ of the consensus peaks,
2. assign weight 0 to control peaks that overlap with consensus peaks,
3. assign weight as the interpolated $D$ of the other control peaks based on their (distance to nearest TSS, GC content).

Do weighted random sampling using the generated weights.

Come back to compare the $D$ of the consensus peaks and sampled matched controls.
