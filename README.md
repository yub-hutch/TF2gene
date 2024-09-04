
## Install

`devtools::install_github('yub-hutch/TF2gene', auth_token = "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")`


## ATAC-seq peaks

- Summarize peaks in terms of distance to TSS & GC content: `summarize_consensus_peaks`

### Real peaks

- Read bed file of consensus peaks: `format_consensus_peaks`

### Control peaks

- Get coordinaltes of null peaks: `get_control_peaks`
- Write control peaks to bed file: `write_control_peak_to_bed`


## Link ATAC-seq peaks to genes

- Calculate distance between peaks & genes: `calc_peak2gene_distance`
