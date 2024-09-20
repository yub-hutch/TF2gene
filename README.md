
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

## 3. Motif info

```
meta_motif = get_cb_motif_info('/fh/fast/sun_w/kenny_zhang/v10nr_clust_public/singletons', ncores = 36)
```

10,249 motif IDs in total, among which 799 contains multiple motifs. The largest motif cluster contains 80,183 motifs.

When calculating Cluster-Buster score for each motif separately, prioritize motif IDs by cluster size. (create_cistarget_motif_databases.py does this iternally when paralleling for multiple motif IDs.)

```
sorted_motifs = meta_motif %>% arrange(desc(num_inner_motif)) %>% pull(motif)

write.table(sorted_motifs, file = 'sorted_motifs.lst', row.names = F, col.names = F, quote = F)
```

## 4. Build a database of Cluster-Buster score on matched control peaks

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

10,249 motifs in total. Do for each motif separately.

```
while IFS= read -r motif; do
  echo "$motif" > "feather/$motif.lst"
done < /fh/fast/sun_w/kenny_zhang/v10nr_clust_public/motifs.lst
```

```
#!/bin/bash

#SBATCH --job-name=CB_control
#SBATCH --output=feather/slurm-%A_%a.out
#SBATCH --error=feather/slurm-%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB

ml Cluster-Buster/0.0-GCC-12.2.0
ml SciPy-bundle/2023.02-gfbf-2022b
ml flatbuffers-python/23.1.4-GCCcore-12.2.0

motif=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sorted_motifs.lst)

echo ${motif} > feather/"${motif}.lst"

python3 /fh/fast/sun_w/kenny_zhang/create_cisTarget_databases/create_cistarget_motif_databases.py \
        -f null_peaks.fa \
        -M /fh/fast/sun_w/kenny_zhang/v10nr_clust_public/singletons \
        -m feather/${motif}.lst \
        -o feather/${motif} \
        -t 1
```

`create_cistarget_motif_databases.py` creates 3 files for each motif:

- motifs_vs_regions.scores.feather (file size ~30M)
- regions_vs_motifs.scores.feather (~ 1G), which is simply the transpose of motifs_vs_regions.scores.feather,
- regions_vs_motifs.rankings.feather (~ 1G), which is the ranking of regions for each motif.

Therefore, keeping only motifs_vs_regions.scores.feather is enough. May supresss writing regions_vs_motifs.scores.feather & regions_vs_motifs.rankings.feather to speed up (by doing some simple edits to the souce code) if more computation is needed in the future.

## 5. Calculate P-value for Cluster-buster score

To sample matched controls of ATAC-seq consensus peaks, generate weights for control peaks in 3 steps:
1. fit (distance to nearest TSS, GC content) joint density $D$ of the consensus peaks,
2. assign weight 0 to control peaks that overlap with consensus peaks,
3. assign weight as the interpolated $D$ of the other control peaks based on their (distance to nearest TSS, GC content).

Do weighted random sampling using the generated weights.

Come back to compare the $D$ of the consensus peaks and sampled matched controls.
