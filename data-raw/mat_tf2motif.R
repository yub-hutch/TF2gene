# # Load data
# tf2motif = fread('/fh/fast/sun_w/kenny_zhang/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl') %>%
#   rename(motif_id = '#motif_id')
#
# dim(tf2motif)
#
#
# # EDA
# types = c('gene is directly annotated',
#           'motif similar to',
#           'gene is orthologous to',
#           'motif is annotated for orthologous gene',
#           'gene is annotated for similar motif')
#
# num_types = sapply(types, function(type) sum(startsWith(tf2motif$description, type)))
#
# sum(num_types) == nrow(tf2motif)
#
# tf2motif$type = sapply(tf2motif$description, function(description) {
#   which_type = which(sapply(types, function(type) startsWith(description, type)))
#   types[which_type]
# })
#
#
# # Calculate TF2motif matrix
# tfs = sort(unique(tf2motif$gene_name))
# motifs = sort(unique(tf2motif$motif_id))
#
# # With direct annotation
# df_direct = tf2motif %>%
#   filter(type == 'gene is directly annotated') %>%
#   select(gene_name, motif_id) %>%
#   unique()
#
# mat_direct = Matrix::sparseMatrix(
#   i = match(df_direct$gene_name, tfs),
#   j = match(df_direct$motif_id, motifs),
#   dims = c(length(tfs), length(motifs)),
#   dimnames = list(tfs, motifs)
# )
#
# # With orthology
# df_orth = tf2motif %>%
#   filter(type %in% c('gene is orthologous to', 'motif is annotated for orthologous gene')) %>%
#   select(gene_name, motif_id) %>%
#   unique()
#
# mat_orth = Matrix::sparseMatrix(
#   i = match(df_orth$gene_name, tfs),
#   j = match(df_orth$motif_id, motifs),
#   dims = c(length(tfs), length(motifs)),
#   dimnames = list(tfs, motifs)
# )
#
#
# # With any annotation
# df_all = tf2motif %>%
#   select(gene_name, motif_id) %>%
#   unique()
#
# mat_simi = Matrix::sparseMatrix(
#   i = match(df_all$gene_name, tfs),
#   j = match(df_all$motif_id, motifs),
#   dims = c(length(tfs), length(motifs)),
#   dimnames = list(tfs, motifs)
# )
#
# # Combine
# mat = pmax(3 * mat_direct, 2 * mat_orth, mat_simi)
