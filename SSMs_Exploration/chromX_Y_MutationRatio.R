
setwd("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/ChromX_Y_MutationRate/")

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

chrX <- integer(length = length(vcf_files)) 
chrY <- integer(length = length(vcf_files))

for (i in 1:length(vcf_files))
{
  #Esto lo tenemos que hacer por el error de la entrada de INFO - temporalmente creamos un archivo sin esa info
  commandRemINFO <- paste0('zgrep -v "INFO=<ID=dbsnp_somatic,Number=.,Type=Flag" ',vcf_files[i], " | gzip -c > tmp.vcf.gz")
  system(commandRemINFO)
  
  #Total number of lines
  all <- system("zcat tmp.vcf.gz | wc -l", intern = TRUE)
  
  #Number of mutations referring to Chromosome X
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr X --recode --stdout | wc -l ", intern=TRUE)
  chrX[i] <- as.integer(all) - as.integer(aux)
  
  #Number of mutations referring to Chromosome Y
  aux <- system("vcftools --gzvcf tmp.vcf.gz --not-chr Y --recode --stdout | wc -l ", intern=TRUE)
  chrY[i] <- as.integer(all) - as.integer(aux)
}

mut_rates_sexChrom <- data.frame(chrX = chrX, chrY = chrY, filename = file_names)
#mut_rates_sexChrom$filename <- as.character(mut_rates_sexChrom$filename) #Forgot to stringAsfactors=FALSE

# En el mismo workspace cargamos info de todos los SSM96s en WGS, además de info de Gender con aliquot_id, etc...
save.image("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/ChromX_Y_MutationRate/Sex_Chrom_MutRates.RData")

# Los objetos t_sample2gender and t_sample2tumortype ya tienen el mismo orden == tot_SSM_WGS
# Lo más fácil será reordenar el df mut_rates_sexChrom como el resto de objetos

# Antes de reordenar vamos a eliminar aquellos que files con los que no trabajamos (blacklisted files)
filenames_sexchrom <- sapply(mut_rates_sexChrom$filename, 
                             function(x) unlist(strsplit(x, split=".consensus.20160830.somatic.snv_mnv.vcf.gz",fixed=TRUE))[1])

new_order <- match(filenames_sexchrom, t_sample2tumortype$tumor_wgs_aliquot_id)
mut_rates_sexChrom <- mut_rates_sexChrom[-(which(is.na(new_order))),]

# Ahora sí podemos buscar el orden adecuado para los 2583 finales
filenames_sexchrom <- sapply(mut_rates_sexChrom$filename, 
                             function(x) unlist(strsplit(x, split=".consensus.20160830.somatic.snv_mnv.vcf.gz",fixed=TRUE))[1])
new_order <- match(t_sample2tumortype$tumor_wgs_aliquot_id, filenames_sexchrom)
mut_rates_sexChrom <- mut_rates_sexChrom[new_order, ]

# Eliminamos columna de nombre de fichero y adjuntamos resto de datos interesantes

mut_rates_sexChrom <- mut_rates_sexChrom[,-3]
mut_rates_sexChrom <- cbind(mut_rates_sexChrom, allchrom = tot_SSM_WGS,
                            Gender = t_sample2gender$donor_sex,
                            Tumor_Type = t_sample2tumortype$histology_abbreviation)
rownames(mut_rates_sexChrom) <- names(tot_SSM_WGS)


# Analysis mutation rate chrX y chrY --------------------------------------

# Expected percentage: based on the proportion of bases in chrX/chrY with respect to all chromosomes (WGS)
# Calculated from WGS_WES_kmer_comparison
perc_chrX_expected <- 5.28
perc_chrY_expected <- 0.897

mut_rates_sexChrom <- cbind(mut_rates_sexChrom, 
                            perc_chrX = (mut_rates_sexChrom$chrX/mut_rates_sexChrom$allchrom)*100, 
                            perc_chrY = (mut_rates_sexChrom$chrY/mut_rates_sexChrom$allchrom)*100)

# Plot percentages per tumor type

# Chromosome X: including the gender
#Let's reorder the tumour types as a function of the percentage median values:
median_values <- ddply(mut_rates_sexChrom,~Tumor_Type,summarise,median = median(perc_chrX,na.rm=TRUE))
order <- median_values$Tumor_Type[order(median_values$median)]
mut_rates_sexChrom$Tumor_Type <- factor(mut_rates_sexChrom$Tumor_Type, levels = order)

boxplot <- ggplot(mut_rates_sexChrom, aes(y=perc_chrX, x=Tumor_Type,fill=Gender)) + 
  geom_boxplot() +
  ggtitle("Percentage of mutations in chrX compared to all genome") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15))) +
  geom_hline(yintercept = perc_chrX_expected, color="red", lty = "dashed") 

pdf("chrX_mutations_perc.pdf",width=12, height=14)
boxplot
dev.off()

# Chromosome Y: only valid for men
#Let's reorder the tumour types as a function of the percentage median values:
median_values <- ddply(mut_rates_sexChrom,~Tumor_Type,summarise,median = median(perc_chrY,na.rm=TRUE))
order <- median_values$Tumor_Type[order(median_values$median)]
mut_rates_sexChrom$Tumor_Type <- factor(mut_rates_sexChrom$Tumor_Type, levels = order)

boxplot <- ggplot(mut_rates_sexChrom, aes(y=perc_chrY, x=Tumor_Type)) + 
  geom_boxplot() +
  ggtitle("Percentage of mutations in chrY compared to all genome") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15))) +
  geom_hline(yintercept = perc_chrY_expected, color="red", lty = "dashed") 

pdf("chrY_mutations_perc.pdf",width=12, height=14)
boxplot
dev.off()
