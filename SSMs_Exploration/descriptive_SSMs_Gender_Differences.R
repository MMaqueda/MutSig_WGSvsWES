# Percentage of mutations WES in WGS per sample ---------------------------

# This information is complementary to the results of the Hypergeometric test
# We are not interested in the 3mer percentages => general percentage based on the
# number of mutations per samples in WES compared to WGS

setwd("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/")
load("Comparison_WGS_WES.RData")

#Compute percentages

tot_SSM_cds <- colSums(mut96_mat_allCDS)
tot_SSM_exon <- colSums(mut96_mat_allExon)
tot_SSM_WGS <- colSums(mut96_mat_allWGS)

perc_WESexon <- (tot_SSM_exon / tot_SSM_WGS)*100
perc_WEScds <- (tot_SSM_cds / tot_SSM_WGS)*100

###########
#If we have to include gender in the samples
load("t_sample2gender.RData")
load("/project/devel/PCAWG/somatic_consensus_aug2016/data/t_sample2tumortype_aliquot.RData")

#Primero hay que reordenar 't_sample2tumortype' tal como está en mut_mat (== names(percWESexon)). Esta
#ordenación se tiene que hacer en base al t_sample.index (el id es único)

sample_index <- sapply(names(perc_WESexon), function(x) unlist(strsplit(x, split=".",fixed=TRUE))[1])
new_order <- match(sample_index, t_sample2tumortype$t_sample_index)
t_sample2tumortype <- t_sample2tumortype[new_order,]

#Una vez reordenado 't_sample2tumortype', ya podemos hacer esto....:
new_order <- match(t_sample2tumortype$tumor_wgs_aliquot_id, t_sample2gender$tumor_wgs_aliquot_id) 
t_sample2gender <- t_sample2gender[new_order,]

toplot_Percentages <- data.frame(WESexon = perc_WESexon, 
                                 WEScds = perc_WEScds, 
                                 Tumour_Type = sapply(names(perc_WESexon), function(x) unlist(strsplit(x, "[.]"))[2]), 
                                 Gender = t_sample2gender$donor_sex,
                                 stringsAsFactors = F)

#Let's reorder the tumour types as a function of the percentage median values:
median_values <- ddply(toplot_Percentages,~Tumour_Type,summarise,median = median(WEScds,na.rm=TRUE))
order <- median_values$Tumour_Type[order(median_values$median)]
toplot_Percentages$Tumour_Type <- factor(toplot_Percentages$Tumour_Type, levels = order)

cols_pcawg <- pcawg.colour.palette(x = levels(toplot_Percentages$Tumour_Type),scheme = "tumour.subtype")
names(cols_pcawg) <- NULL

give.n <- function(x){
  return(c(y = min(x)*0.85, label = length(x))) 
}

require("ggpubr")

boxplot <- ggplot(toplot_Percentages, aes(y=WEScds, x=Tumour_Type,fill=Gender)) + 
  geom_boxplot() +
  #scale_fill_manual(values= cols_pcawg) +
  ggtitle("Percentage number of mutations (SSM96 counts) \n CDS REGION compared to WGS") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15)))+
  geom_hline(yintercept = 1.23, color="red") +
  stat_summary(fun.data = give.n, geom = "text",
               position = position_dodge(width = 0.75),
               size=4, angle=90) 

pdf("Mutations_Percentage_wGENDER_CDS.pdf",width=12, height=14)
boxplot
dev.off()

save.image("Comparison_WGS_WES.RData")

boxplot <- ggplot(toplot_Percentages, aes(y=WEScds, x=Gender)) +
  facet_wrap( ~ Tumour_Type) +
  geom_boxplot() +
  #scale_fill_manual(values= cols_pcawg) +
  ggtitle("Percentage number of mutations (SSM96 counts) \n CDS REGION compared to WGS") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15)))+
  geom_hline(yintercept = 1.23, color="red") +
  stat_summary(fun.data = give.n, geom = "text",
               position = position_dodge(width = 0.75),
               size=4, angle=90) + 
  stat_compare_means()

pdf("Mutations_Percentage_wGENDER_CDSXXX.pdf",width=12, height=14)
boxplot
dev.off()

####################

toplot_totSSMs <- data.frame(WGS = tot_SSM_WGS, 
                             WEScds = tot_SSM_cds,
                             WESexon = tot_SSM_exon,
                             Tumour_Type = sapply(names(tot_SSM_exon), function(x) unlist(strsplit(x, "[.]"))[2]), 
                             Gender = t_sample2gender$donor_sex,
                             stringsAsFactors = F)

#Let's reorder the tumour types as a function of the percentage median values:
median_values <- ddply(toplot_totSSMs,~Tumour_Type,summarise,median = median(WEScds,na.rm=TRUE))
order <- median_values$Tumour_Type[order(median_values$median)]
toplot_totSSMs$Tumour_Type <- factor(toplot_totSSMs$Tumour_Type, levels = order)

boxplot <- ggplot(toplot_totSSMs, aes(y=WEScds, x=Tumour_Type,fill=Gender)) + 
  geom_boxplot() + scale_y_continuous(trans = "log10") +
  ggtitle("Total number of mutations (SSM96 counts) \n WES - CDS data") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15))) 

pdf("Total_SSMs_CDSwGender.pdf",width=12, height=14)
boxplot
dev.off()


boxplot <- ggplot(toplot_totSSMs, aes(y=WEScds, x=Gender, fill=Gender)) +
  facet_wrap( ~ Tumour_Type, scales="free") +
  geom_boxplot() + scale_y_continuous(trans = "log10") +
  ggtitle("Total number of mutations (SSM96 counts) \n WES CDS data") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15))) +
  stat_summary(fun.data = give.n, geom = "text",
               position = position_dodge(width = 1),
               size=4, angle=90) + 
  stat_compare_means(label="p.signif") 

pdf("xx.pdf",width=12, height=14)
boxplot
dev.off()

# Now we should correct the number of chrX mutations for male samples
# We'll only do that for the WGS data since we only have the mutations in WGS (not in WES) 

load("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/ChromX_Y_MutationRate/Sex_Chrom_MutRates.RData")

tosum <- as.integer()
for(i in 1:dim(mut_rates_sexChrom)[1])
{
  if(mut_rates_sexChrom$Gender[i] == "female") {tosum[i] <- 0}
  else if(mut_rates_sexChrom$Gender[i] == "male") {tosum[i] <- mut_rates_sexChrom$chrX[i]}
}

toplot_totSSMs$WGS_corrGender <- toplot_totSSMs$WGS + tosum

boxplot <- ggplot(toplot_totSSMs, aes(y=WGS_corrGender, x=Gender, fill=Gender)) +
  facet_wrap( ~ Tumour_Type) +
  geom_boxplot() + scale_y_continuous(trans = "log10") +
  ggtitle("Comparison SSMs Female vs Male in WGS data \ chrX in Male is being corrected (x2)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15))) +
  stat_summary(fun.data = give.n, geom = "text",
               position = position_dodge(width = 1),
               size=4, angle=90) + 
  stat_compare_means(label="p.signif") 

pdf("WGS_Comparison_chrXMale_corrected.pdf",width=12, height=14)
boxplot
dev.off()  

  


