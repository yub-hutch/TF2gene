# control_region = generate_control_regions(meta_gene = grch38, chr_lens = grch38_chr_lens, region_length = 500, radius = 500e3)
# write_control_regions_to_bed(control_region, fbed = 'control_regions.bed')
#
#
# ml BEDTools/2.30.0-GCC-12.2.0
#
# bedtools getfasta \
#   -fi ~/.nextflow/assets/nf-core/grch38_genome/genome.fa \
#   -bed output/control_regions/control_regions.bed \
#   -fo output/control_regions/control_regions.fa
#
#
# meta_control_region = summarize_features_of_regions(fasta = 'control_regions.fa', ncores = 36, meta_gene = TF2gene::grch38)
