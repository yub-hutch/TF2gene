# Downloaded from https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl.

# tf2motif = fread('/fh/fast/sun_w/kenny_zhang/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl') %>%
#   rename(motif_id = '#motif_id') %>%
#   as_tibble()
#
# tf2motif$type = sapply(tf2motif$description, function(description) {
#   which_type = which(sapply(types, function(type) startsWith(description, type)))
#   types[which_type]
# })
