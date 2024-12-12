
## 1. Install

`devtools::install_github('yub-hutch/TF2gene')`


## 2. Usage

```
result = TF2gene(
  feather_file,
  fasta_consensus_peaks,
  dir_null_cbscore,
  peak2gene,
  dir_save,
  resume,
  num_control_per_peak = 1000,
  rules = c(0.05, 0.01, 0.25)
)
```

where
- `feather_file` is the feather file storing the Cluster-Buster scores on consensus peaks,
- `fasta_consensus_peaks` is the sequence file of the consenesus peaks,
- `dir_null_cbscore` is the directory storing Cluster-Buster scrores on genome-wide control regions. Use `/fh/fast/sun_w/yub/grn/rpe1/annotation/output/control_regions/feather` for now.
- `peak2gene` is the sparse binary matrix representing the connection between consensus peaks and genes,
- `dir_save` is the directory to save intermediate outputs, which will be loaded when the unfinished task is resumed,
- `resume` is logical indicating whether to resume the unfinished task,
- `num_control_per_peak` is the number of matched control regions to use for each consensus peak,
- `rules` is used for motif selection. (1) Proportion(P-value < 0.05) > 5%, (2) Proportion(P-value < 0.01) > 1%, (3) Minimum Q-value < 0.25. Tune `rules` to get a reasonable number of active motifs.


To get `peak2gene` by distance, use `get_peak2gene`.


## 3. Notes

### Motifs

```
meta_motif = get_cb_motif_info('/fh/fast/sun_w/kenny_zhang/v10nr_clust_public/singletons', ncores = 36)
```

10,249 motif IDs in total, among which 799 contains multiple motifs. The largest motif cluster contains 80,183 motifs.

When calculating Cluster-Buster score for each motif separately, prioritize motif IDs by cluster size. (create_cistarget_motif_databases.py does this iternally when paralleling for multiple motif IDs.)

```
sorted_motifs = meta_motif %>% arrange(desc(num_inner_motif)) %>% pull(motif)

write.table(sorted_motifs, file = 'sorted_motifs.lst', row.names = F, col.names = F, quote = F)
```

### Build a database of Cluster-Buster score on matched control peaks

Run Cluster-Buster on control peaks:

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
