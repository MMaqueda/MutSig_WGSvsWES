setwd("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/Gral_Chrom_Rate/")

#Los archivos de los que queremos saber el número de mutaciones: en un principio lo hacemos solo para WGS
VCFfilesPath = "/project/devel/PCAWG/somatic_consensus_aug2016/data/filtered_calls/whitelist_12oct/ssm"
#VCFfilesPath = "/project/devel/PCAWG/mmaqueda/WES_VCFfiles/VCF_Filtered/CDS"
#VCFfilesPath = "/project/devel/PCAWG/mmaqueda/WES_VCFfiles/VCF_Filtered/EXON"

vcf_files <- list.files(path=VCFfilesPath,pattern = ".vcf.gz", full.names = TRUE) #Including also .tbi
#Excluding .tbi files (valid for WGS ONLY)
vcf_files <- vcf_files[-(grep(pattern=".tbi", x= vcf_files, fixed=TRUE))]

#Los archivos generados una vez filtrados:
file_names <- list.files(path=VCFfilesPath,pattern = ".vcf.gz", full.names = FALSE)
file_names <- file_names[-(grep(pattern= ".tbi", x=file_names, fixed=TRUE))]

chr1 <- integer(length = length(vcf_files)) 
chr2 <- integer(length = length(vcf_files))
chr3 <- integer(length = length(vcf_files))
chr4 <- integer(length = length(vcf_files))
chr5 <- integer(length = length(vcf_files))
chr6 <- integer(length = length(vcf_files))
chr7 <- integer(length = length(vcf_files))
chr8 <- integer(length = length(vcf_files))
chr9 <- integer(length = length(vcf_files))
chr10 <- integer(length = length(vcf_files))
chr11 <- integer(length = length(vcf_files))
chr12 <- integer(length = length(vcf_files))
chr13 <- integer(length = length(vcf_files))
chr14 <- integer(length = length(vcf_files))
chr15 <- integer(length = length(vcf_files))
chr16 <- integer(length = length(vcf_files))
chr17 <- integer(length = length(vcf_files))
chr18 <- integer(length = length(vcf_files))
chr19 <- integer(length = length(vcf_files))
chr20 <- integer(length = length(vcf_files))
chr21 <- integer(length = length(vcf_files))
chr22 <- integer(length = length(vcf_files))

for (i in 1:length(vcf_files))
{
  #Esto lo tenemos que hacer por el error de la entrada de INFO - temporalmente creamos un archivo sin esa info
  commandRemINFO <- paste0('zgrep -v "INFO=<ID=dbsnp_somatic,Number=.,Type=Flag" ',vcf_files[i], " | gzip -c > tmp.vcf.gz")
  system(commandRemINFO)
  
  #Total number of lines
  all <- system("zcat tmp.vcf.gz | wc -l", intern = TRUE)
  
  #Number of mutations referring to Chromosome 1
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 1 --recode --stdout | wc -l ", intern=TRUE)
  chr1[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 2
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 2 --recode --stdout | wc -l ", intern=TRUE)
  chr2[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 3
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 3 --recode --stdout | wc -l ", intern=TRUE)
  chr3[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 4
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 4 --recode --stdout | wc -l ", intern=TRUE)
  chr4[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 5
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 5 --recode --stdout | wc -l ", intern=TRUE)
  chr5[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 6
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 6 --recode --stdout | wc -l ", intern=TRUE)
  chr6[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 7
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 7 --recode --stdout | wc -l ", intern=TRUE)
  chr7[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 8
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 8 --recode --stdout | wc -l ", intern=TRUE)
  chr8[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 9
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 9 --recode --stdout | wc -l ", intern=TRUE)
  chr9[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 10
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 10 --recode --stdout | wc -l ", intern=TRUE)
  chr10[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 11
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 11 --recode --stdout | wc -l ", intern=TRUE)
  chr11[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 12
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 12 --recode --stdout | wc -l ", intern=TRUE)
  chr12[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 13
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 13 --recode --stdout | wc -l ", intern=TRUE)
  chr13[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 14
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 14 --recode --stdout | wc -l ", intern=TRUE)
  chr14[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 15
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 15 --recode --stdout | wc -l ", intern=TRUE)
  chr15[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 16
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 16 --recode --stdout | wc -l ", intern=TRUE)
  chr16[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 17
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 17 --recode --stdout | wc -l ", intern=TRUE)
  chr17[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 18
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 18 --recode --stdout | wc -l ", intern=TRUE)
  chr18[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 19
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 19 --recode --stdout | wc -l ", intern=TRUE)
  chr19[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 20
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 20 --recode --stdout | wc -l ", intern=TRUE)
  chr20[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 21
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 21 --recode --stdout | wc -l ", intern=TRUE)
  chr21[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome 22
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr 22 --recode --stdout | wc -l ", intern=TRUE)
  chr22[i] <- as.integer(all) - as.integer(aux)
  
}

mut_rates_autosChrom <- data.frame(chr1 = chr1, chr2 = chr2, chr3 = chr3, chr4 = chr4,
                                   chr5 = chr5, chr6 = chr6, chr7 = chr7, chr8 = chr8,
                                   chr9 = chr9, chr10 = chr10, chr11 = chr11,
                                   chr12 = chr12, chr13 = chr13, chr14 = chr14,
                                   chr15 = chr15, chr16 = chr16, chr17 = chr17,
                                   chr18 = chr18, chr19= chr19, chr20 = chr20,
                                   chr21 = chr21, chr22 = chr22,
                                 filename = file_names, stringsAsFactors = FALSE)

# En el mismo workspace cargamos info de todos los SSM96s en WGS, además de info de Gender con aliquot_id, etc...
save.image("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/Gral_Chrom_Rate/Autosomal_Chrom_MutRates.RData")

# Los objetos t_sample2gender and t_sample2tumortype ya tienen el mismo orden == tot_SSM_WGS
# Lo más fácil será reordenar el df mut_rates_sexChrom como el resto de objetos

# Antes de reordenar vamos a eliminar aquellos que files con los que no trabajamos (blacklisted files)
# filenames_sexchrom <- sapply(mut_rates_sexChrom$filename, 
#                              function(x) unlist(strsplit(x, split=".consensus.20160830.somatic.snv_mnv.vcf.gz",fixed=TRUE))[1])
# 
# new_order <- match(filenames_sexchrom, t_sample2tumortype$tumor_wgs_aliquot_id)
# mut_rates_sexChrom <- mut_rates_sexChrom[-(which(is.na(new_order))),]
# 
# # Ahora sí podemos buscar el orden adecuado para los 2583 finales
# filenames_sexchrom <- sapply(mut_rates_sexChrom$filename, 
#                              function(x) unlist(strsplit(x, split=".consensus.20160830.somatic.snv_mnv.vcf.gz",fixed=TRUE))[1])
# new_order <- match(t_sample2tumortype$tumor_wgs_aliquot_id, filenames_sexchrom)
# mut_rates_sexChrom <- mut_rates_sexChrom[new_order, ]
# 
# # Eliminamos columna de nombre de fichero y adjuntamos resto de datos interesantes
# 
# mut_rates_sexChrom <- mut_rates_sexChrom[,-3]
# mut_rates_sexChrom <- cbind(mut_rates_sexChrom, allchrom = tot_SSM_WGS,
#                             Gender = t_sample2gender$donor_sex,
#                             Tumor_Type = t_sample2tumortype$histology_abbreviation)
# rownames(mut_rates_sexChrom) <- names(tot_SSM_WGS)
# 
# 
# # Analysis mutation rate chrX y chrY --------------------------------------
# 
# # Expected percentage: based on the proportion of bases in chrX/chrY with respect to all chromosomes (WGS)
# # Calculated from WGS_WES_kmer_comparison
# perc_chrX_expected <- 5.28
# perc_chrY_expected <- 0.897
# 
# mut_rates_sexChrom <- cbind(mut_rates_sexChrom, 
#                             perc_chrX = (mut_rates_sexChrom$chrX/mut_rates_sexChrom$allchrom)*100, 
#                             perc_chrY = (mut_rates_sexChrom$chrY/mut_rates_sexChrom$allchrom)*100)
# 
# # Plot percentages per tumor type
# 
# # Chromosome X: including the gender
# #Let's reorder the tumour types as a function of the percentage median values:
# median_values <- ddply(mut_rates_sexChrom,~Tumor_Type,summarise,median = median(perc_chrX,na.rm=TRUE))
# order <- median_values$Tumor_Type[order(median_values$median)]
# mut_rates_sexChrom$Tumor_Type <- factor(mut_rates_sexChrom$Tumor_Type, levels = order)
# 
# boxplot <- ggplot(mut_rates_sexChrom, aes(y=perc_chrX, x=Tumor_Type,fill=Gender)) + 
#   geom_boxplot() +
#   ggtitle("Percentage of mutations in chrX compared to all genome") +
#   theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
#                      axis.text.y = element_text(size=14),
#                      plot.title = element_text(face = "bold", size = (15))) +
#   geom_hline(yintercept = perc_chrX_expected, color="red", lty = "dashed") 
# 
# pdf("chrX_mutations_perc.pdf",width=12, height=14)
# boxplot
# dev.off()
# 
# # Chromosome Y: only valid for men
# #Let's reorder the tumour types as a function of the percentage median values:
# median_values <- ddply(mut_rates_sexChrom,~Tumor_Type,summarise,median = median(perc_chrY,na.rm=TRUE))
# order <- median_values$Tumor_Type[order(median_values$median)]
# mut_rates_sexChrom$Tumor_Type <- factor(mut_rates_sexChrom$Tumor_Type, levels = order)
# 
# boxplot <- ggplot(mut_rates_sexChrom, aes(y=perc_chrY, x=Tumor_Type)) + 
#   geom_boxplot() +
#   ggtitle("Percentage of mutations in chrY compared to all genome") +
#   theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
#                      axis.text.y = element_text(size=14),
#                      plot.title = element_text(face = "bold", size = (15))) +
#   geom_hline(yintercept = perc_chrY_expected, color="red", lty = "dashed") 
# 
# pdf("chrY_mutations_perc.pdf",width=12, height=14)
# boxplot
# dev.off()
