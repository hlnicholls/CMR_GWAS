library(CMplot)
library(data.table)
library(dplyr)

lv_phenotypes <- c("LVEDV", "LVESV", "LVSV", "LVEF", "LVM", "LVMVR", "LVGFI", "LVMCF", "LVEDV_BSA", "LVESV_BSA", "LVSV_BSA", "LVM_BSA")
rv_phenotypes <- c("RVEDV", "RVESV", "RVSV", "RVEF", "RV_LV_ratio", "RVEDV_BSA", "RVESV_BSA", "RVSV_BSA")

lv_df <- data.frame()

for(pheno in lv_phenotypes) {
  file_path <- paste0('/Input/', pheno, '_GWAS_38_37.txt')
  temp_df <- fread(file_path, select = c("CHROM", "GENPOS", "SNP", "P"))
  setnames(temp_df, "P", pheno)
  if(nrow(lv_df) == 0) {
    lv_df <- temp_df
  } else {
    lv_df <- merge(lv_df, temp_df, all = TRUE)
  }
}

setnames(lv_df, c('SNP', 'CHROM', 'GENPOS'), c('SNP', 'Chromosome', 'Position'))
lv_df <- lv_df[, c("SNP", "Chromosome", "Position", "LVEDV", "LVESV", "LVSV", "LVEF", "LVM", "LVMVR", "LVGFI", "LVMCF", "LVEDV_BSA", "LVESV_BSA", "LVSV_BSA", "LVM_BSA")]

lv_df <- lv_df[!duplicated(lv_df$SNP), ]

setwd('/Output/Plots')

# Checking possible colours
#library(pals)
#stepped_colors <- stepped(24)
#stepped

# Define a unique color for each LV trait
lv_colors <- matrix(c(
  "#990F26", "#B33E52", "#CC7A88", "#E6B8BF","#0F8299", "#3E9FB3", "#7ABECC", "#B8DEE6",
  "#3D0F99", "#653EB3", "#967ACC", "#C7B8E6"
), nrow = length(lv_phenotypes), byrow = TRUE)


CMplot(lv_df,
       type = "p",
       plot.type = "c",
       chr.labels = paste("Chr", c(1:22), sep = ""),
       col = lv_colors,
       r = 0.4,
       cir.axis = TRUE,
       outward = FALSE,
       cir.axis.col = "black",
       cir.chr = TRUE, 
       cir.chr.h = 1.3,
       chr.den.col = "black",
       file = "jpg",
       file.name = "01AUG24circ_all_lv_traits_multicolour",
       dpi = 300,
       file.output = TRUE,
       verbose = TRUE,
       width = 15,
       height = 15)

rv_df <- data.frame()

for(pheno in rv_phenotypes) {
  file_path <- paste0('/Input/', pheno, '_GWAS_38_37.txt')
  temp_df <- fread(file_path, select = c("CHROM", "GENPOS", "SNP", "P"))
  setnames(temp_df, "P", pheno)
  if(nrow(rv_df) == 0) {
    rv_df <- temp_df
  } else {
    rv_df <- merge(rv_df, temp_df, all = TRUE)
  }
}

setnames(rv_df, c('SNP', 'CHROM', 'GENPOS'), c('SNP', 'Chromosome', 'Position'))
rv_df <- rv_df[, c("SNP", "Chromosome", "Position", "RVEDV", "RVESV", "RVSV", "RVEF", "RV_LV_ratio", "RVEDV_BSA", "RVESV_BSA", "RVSV_BSA")]

rv_df <- rv_df[!duplicated(rv_df$SNP), ]

# Define a unique color for each RV trait
rv_colors <- matrix(c(
  "#99600F", "#CCAA7A", "#E6D2B8", "#54990F",  "#A3CC7A", "#CFE6B8","#333333", "#CCCCCC"
), nrow = length(rv_phenotypes), byrow = TRUE)

CMplot(rv_df,
       type = "p",
       plot.type = "c",
       chr.labels = paste("Chr", c(1:22), sep = ""),
       col = rv_colors,
       r = 0.4,
       cir.axis = TRUE,
       outward = FALSE,
       cir.axis.col = "black",
       cir.chr = TRUE,
       cir.chr.h = 1.3,
       chr.den.col = "black",
       file = "jpg",
       file.name = "01AUG24circ_all_rv_traits_multicolour",
       dpi = 300,
       file.output = TRUE,
       verbose = TRUE,
       width = 15,
       height = 15)
