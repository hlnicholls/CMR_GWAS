library(topr)
library(data.table)
library(dplyr)

# Repeat code for LV phenotypes
phenotypes <- c("RVEDV", "RVESV", "RVSV", "RVEF", "RV_LV_ratio", "RVEDV_BSA", "RVESV_BSA", "RVSV_BSA")

for (phenotype in phenotypes) {

  file_path <- paste0('/Input/', phenotype, '_assoc_regenie_allchr.txt')
  print('reading in file')
  gwas <- fread(file_path, select = c('CHROM', 'GENPOS', 'ALLELE0', 'ALLELE1', 'p', 'A1FREQ', 'BETA', 'SE'))

  colnames(gwas)[c(2, 3, 4, 5, 6)] <- c('POS', 'REF', 'ALT', 'P', 'AF')

  gwas_filtered <- gwas[gwas$P < 5e-8, ]
  print('get genes')
  gwas_annotated <- annotate_with_nearest_gene(gwas_filtered)
  gwas[gwas$P < 5e-8, "Gene_Symbol"] <- gwas_annotated$Gene_Symbol
  print('plotting...')
  png_filename <- paste0("/Output/topr_manhattan/manhattan38_", phenotype, "_annotated.png")
  png(png_filename, width = 3000, height = 1800, res = 300)
  plot(manhattan(gwas, annotate = 5e-9, title=phenotype))
  dev.off()
}
