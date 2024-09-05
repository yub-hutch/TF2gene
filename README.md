
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

5,807,331 control peaks in total.

### 2. Summarize control peaks

```
meta_null_peak = summarize_consensus_peaks(fasta = 'null_peaks.fa', ncores = 36)

saveRDS(meta_null_peak, file = 'meta_null_peak.rds')
```

```
library(ggpubr)

ggplot(meta_null_peak, aes(dist, gc_content)) +
  geom_hex(bins = 100) +
  labs(x = 'Distance to nearest TSS', y = 'GC content') +
  stat_cor(method = 'spearman', label.x.npc = 'center') +
  guides(fill = 'none') +
  theme_pubr()
```

![image](https://github.com/user-attachments/assets/fc180974-f8f8-4e0e-960a-002805a9a711)

### 3. Run Cluster-Buster on control peaks

Split into 10 batches.

```
total_sequences=$(grep -c "^>" null_peaks.fa)
sequences_per_file=$((total_sequences / 10))

awk -v seqs_per_file="$sequences_per_file" -v total_sequences="$total_sequences" '
BEGIN {
    file_num = 1;
    seq_count = 0;
    out = "null_peaks_batch" file_num ".fa";
    print "Processing file: " out;
}
/^>/ {
    if (seq_count % seqs_per_file == 0 && seq_count != 0 && file_num < 10) {
        close(out);
        file_num++;
        out = "null_peaks_batch" file_num ".fa";
        print "Processing file: " out;
    }
    seq_count++;
    if (seq_count % 10000 == 0) {
            printf "Processed %d of %d sequences (%.2f%%)\n", seq_count, total_sequences, (seq_count / total_sequences) * 100;
        }
    }
{
    print > out
}
END {
    close(out)
    print "Completed processing all sequences.";
}' null_peaks.fa
```

```
#!/bin/bash
#SBATCH --job-name=CB_on_controls
#SBATCH --output=CB_on_controls_output_%A_%a.txt
#SBATCH --error=CB_on_controls_error_%A_%a.txt
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=128GB
#SBATCH --array=1-10

ml Cluster-Buster/0.0-GCC-12.2.0
ml SciPy-bundle/2023.02-gfbf-2022b
ml flatbuffers-python/23.1.4-GCCcore-12.2.0

input="null_peaks_batch${SLURM_ARRAY_TASK_ID}.fa"
out_prefix="null_peaks_batch${SLURM_ARRAY_TASK_ID}"

python3 /fh/fast/sun_w/kenny_zhang/create_cisTarget_databases/create_cistarget_motif_databases.py \
    -f $input \
    -M /fh/fast/sun_w/kenny_zhang/v10nr_clust_public/singletons \
    -m /fh/fast/sun_w/kenny_zhang/v10nr_clust_public/motifs.lst \
    -o $out_prefix \
    -t 36
```

