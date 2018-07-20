
source("/project/devel/PCAWG/mmaqueda/Rscripts/pcawg.colour.palette.R")
source("/project/devel/PCAWG/mmaqueda/Rscripts/normalize_ssm96mat.R")
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(SomaticSignatures))

setwd("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/")
#save.image("./Comparison_WGS_WES.RData")
load("./Comparison_WGS_WES.RData")

#Descriptive Statistics about the mutational profile of the samples when considering CDS or Exon
#First, if we wanna retrieve some statistics of all tumour types, we have to combine all txt files


# Files manipulation ------------------------------------------------------

#List all the files
txt_files <- list.files(path="./CDS_based",pattern = ".txt", full.names = TRUE)
# Read the files in
txt_files_df <- lapply(txt_files, function(x) {read.table(file = x, header = T, check.names=F)})
# Combine them
mut96_mat_allCDS <- do.call("cbind", lapply(txt_files_df, as.data.frame)) 

#We can save the df in a txt file in case it is needed
write.table(mut96_mat_allCDS,file = "./mut96_mat_allCDS.txt",sep="\t")

#NOTE: the same process has been repeated for Exon based files

#INPUT: a mutation matrix (SSM96) to be read from a text file in case it is not in the workspace
#mut96_mat <- as.matrix(read.table("./mut96_mat_allCDS.txt",header=TRUE,check.names=FALSE))

#####1 Summary table 1
#Refer to descriptive_SSMs_WGS.R for summary table (tab1): Summary_SSMs_pertumour_CDS
#This includes the ....Preparing the matrix information...


#####2  Figure dotplot with number of counts
#Refer to Figure 8 from descriptive_SSMs_WGS.R to obtain this boxplot



# Hypergeometric test #muts in WES region ---------------------------------

#####Compute the probability of selecting SSM muts in CDS region is lower than random selection from the WGS
#Based on an hypergeometric test.

#First we construct a data frame with all the relevant information for this test
load("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Comparison_Hypergeometric.RData")

#Samples have to follow the same order in the different mutational matrix (CDS, Exon and WGS)
order <- match(colnames(mut96_mat_allWGS),colnames(mut96_mat_allCDS))  
mut96_mat_allCDS <- mut96_mat_allCDS[,order]
order <- match(colnames(mut96_mat_allWGS),colnames(mut96_mat_allExon))  
mut96_mat_allExon <- mut96_mat_allExon[,order]

#And now the df

posWGS <- 2861327131
percentageCDS <- 1.23/100
percentageExon <- 4.26/100
numberSamples <- dim(mut96_mat_allCDS)[2]

Data_hyper <- data.frame(x1SSMsCDS = colSums(mut96_mat_allCDS),
                         x2SSMsExon = colSums(mut96_mat_allExon),
                         mSSMsWGS = colSums(mut96_mat_allWGS),
                         nWGS_woSSMsWGS = posWGS - colSums(mut96_mat_allWGS),
                         k1CDS = rep(posWGS * percentageCDS,numberSamples), 
                         k2Exon = rep(posWGS * percentageExon,numberSamples))

#Let's compute the p-val
#H0: prob of selecting this number of mutations (i.e. CDS region) is no lower that random selection from the WGS
#H1: prob of selecting this number of mutations (i.e. CDS region) is lower that random selection from the WGS



#CDS
adjpval_CDS <- apply(Data_hyper, 1, function(x) p.adjust(phyper(x[1],x[3],x[4],x[5]),"fdr"))  #Depletion
adjpval_CDS_enrich <- apply(Data_hyper, 1, function(x) p.adjust(phyper(x[1],x[3],x[4],x[5],
                                                                       lower.tail = FALSE),"fdr")) #Enrichment

#Exon
adjpval_Exon <- apply(Data_hyper, 1, function(x) p.adjust(phyper(x[2],x[3],x[4],x[6]),"fdr")) #Depletion
adjpval_Exon_enrich <- apply(Data_hyper, 1, function(x) p.adjust(phyper(x[2],x[3],x[4],x[6],
                                                                        lower.tail = FALSE),"fdr")) #Enrichment

#Let's plot the results
#Gathering the needed data
toplot_padj <- data.frame(AdjPvalCDS = adjpval_CDS, 
                     AdjPvalExon = adjpval_Exon, 
                     Tumour_Type = sapply(names(adjpval_CDS_enrich), function(x) unlist(strsplit(x, "[.]"))[2]), 
                     stringsAsFactors = F)

toplot_padj_enrich <- data.frame(AdjPvalCDS = adjpval_CDS_enrich, 
                          AdjPvalExon = adjpval_Exon_enrich, 
                          Tumour_Type = sapply(names(adjpval_CDS_enrich), function(x) unlist(strsplit(x, "[.]"))[2]), 
                          stringsAsFactors = F)

#Let's reorder the tumour types as a function of the adj pval median values:
median_values <- ddply(toplot_padj_enrich,~Tumour_Type,summarise,median = median(AdjPvalCDS,na.rm=TRUE))
order <- median_values$Tumour_Type[order(median_values$median)]
toplot_padj_enrich$Tumour_Type <- factor(toplot_padj$Tumour_Type, levels = order)

cols_pcawg <- pcawg.colour.palette(x = levels(toplot_padj_enrich$Tumour_Type),scheme = "tumour.subtype")
names(cols_pcawg) <- NULL

boxplot <- ggplot(toplot_padj_enrich, aes(y=AdjPvalCDS, x=Tumour_Type, fill= Tumour_Type)) + 
  geom_boxplot() +
  scale_fill_manual(values= cols_pcawg) +
  ggtitle("ENRICHMENT Hypergeometric test - SSMs counts \n CDS REGION") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15)))+
  geom_hline(yintercept = 0.05, color="red")


pdf("Enrichment_Hypergeometric_CDS.pdf",width=12, height=14)
boxplot
dev.off()

# Cosine similarity  mut profile WGS vs WES -------------------------------

#Calculate the cosine similiraty betwen the WGS and WES (CDS or Exon) profile per sample. Results will be plotted
#per tumour type (boxplot of the cosine similarity and )


# No correction considered ------------------------------------------------

#########This is for counts! (or freq is the same)
#For CDS
# calculate all pairwise cosine similarities
cos_sim_WGS_WEScds <- cos_sim_matrix(mut96_mat_allWGS, mut96_mat_allCDS)
# extract cosine similarities per sample between WGS and WES region
cos_sim_WGS_WEScds <- as.data.frame(diag(cos_sim_WGS_WEScds))

#For Exon
# calculate all pairwise cosine similarities
cos_sim_WGS_WESexon <- cos_sim_matrix(mut96_mat_allWGS, mut96_mat_allExon)
# extract cosine similarities per sample between WGS and WES region
cos_sim_WGS_WESexon <- as.data.frame(diag(cos_sim_WGS_WESexon))

#Plot results
colnames(cos_sim_WGS_WEScds) = "cos_sim"
cos_sim_WGS_WEScds$Tumour_Type = sapply(rownames(cos_sim_WGS_WEScds), function(x) unlist(strsplit(x, "[.]"))[2])

colnames(cos_sim_WGS_WESexon) = "cos_sim"
cos_sim_WGS_WESexon$Tumour_Type = sapply(rownames(cos_sim_WGS_WESexon), function(x) unlist(strsplit(x, "[.]"))[2])

#Let's re-order the cos_sim values based on the median cosine similarity of each tumour type
median_values <- ddply(cos_sim_WGS_WEScds,~Tumour_Type,summarise,median = median(cos_sim,na.rm=TRUE))
order <- median_values$Tumour_Type[order(median_values$median)]
cos_sim_WGS_WEScds$Tumour_Type <- factor(cos_sim_WGS_WEScds$Tumour_Type, levels = order)

cols_pcawg <- pcawg.colour.palette(x = levels(cos_sim_WGS_WEScds$Tumour_Type),scheme = "tumour.subtype")
names(cols_pcawg) <- NULL

median_values <- ddply(cos_sim_WGS_WESexon,~Tumour_Type,summarise,median = median(cos_sim,na.rm=TRUE))
order <- median_values$Tumour_Type[order(median_values$median)]
cos_sim_WGS_WESexon$Tumour_Type <- factor(cos_sim_WGS_WESexon$Tumour_Type, levels = order)

cols_pcawg <- pcawg.colour.palette(x = cos_sim_WGS_WESexon$Tumour_Type,scheme = "tumour.subtype")
names(cols_pcawg) <- NULL

boxplot <- ggplot(cos_sim_WGS_WESexon, aes(y=cos_sim, x=Tumour_Type, fill=Tumour_Type)) + 
  geom_boxplot() +
  scale_fill_manual(values = cols_pcawg) +
  ggtitle("Cosine Similarity - SSM96 counts profile \n EXON REGION") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15)))+
  geom_hline(yintercept = 0.9, color="green") +
  geom_hline(yintercept = 0.5, color="red")


pdf("CosSimWGSExoncounts.pdf",width=12, height=14)
boxplot
dev.off()

##############
#This is for frequencies: The result will be exactly the same as in counts since cos sim is avaluating shape not magnitude
#In any case we need the normalized profiles for correction
#########

mut96_mat_allCDS_norm <- normalize_ssm96mat(Path_mut96_mat = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_allCDS.txt", 
                   resultsPath = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/")

mut96_mat_allExon_norm <- normalize_ssm96mat(Path_mut96_mat = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/mut96_mat_allExon.txt", 
                                            resultsPath = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/")

mut96_mat_allWGS_norm <- normalize_ssm96mat(Path_mut96_mat = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/mut96_mat_allWGS.txt", 
                                             resultsPath = "/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/")



#  kmer Correction considered ---------------------------------------------------
###########
# And now we finally compared the mutational profile once WES is corrected to WGS frequencies
# For doing this we will use a function from SomaticSignatures package
require(SomaticSignatures)

#In order to obtain the correction factors, first we need to handle the k3mer frequencies
reduce_k3mer_occur <- function(freq_3mer_WXX)
{
  aux <- colSums(freq_3mer_WXX)  #Originally counts are distributed per chromosome or per region
  #And now the percentage
  k3mer_WXX <- sapply(aux, function(x) x/sum(aux))
  return(k3mer_WXX)
}

k3mer_WGS <- reduce_k3mer_occur(freq_3mer_WGS)
k3mer_CDS <- reduce_k3mer_occur(freq_3mer_CDS)
k3mer_Exon <- reduce_k3mer_occur(freq_3mer_Exon)

#And now we can compute the correction factors:
corr_factorCDS <- k3mer_WGS / k3mer_CDS
corr_factorExon <- k3mer_WGS / k3mer_Exon

#In order to use the correction function we need to renamed the rownames as in the SomSign package

subs <-  c('CA','CG','CT','TA','TC','TG')
context <- c('A.A','A.C','A.G','A.T','C.A','C.C','C.G','C.T',
             'G.A','G.C','G.G','G.T','T.A','T.C','T.G','T.T')

ssm96_somsign <- paste(rep(subs,each=16), rep(context,times=4), sep=" ")
ssm96_std_rownames <- rownames(mut96_mat_allCDS_norm)

rownames(mut96_mat_allCDS_norm) <- ssm96_somsign
rownames(mut96_mat_allExon_norm) <- ssm96_somsign

mut96_mat_allCDS_corr <- normalizeMotifs(mut96_mat_allCDS_norm, corr_factorCDS)
mut96_mat_allExon_corr <- normalizeMotifs(mut96_mat_allExon_norm, corr_factorExon)

#And now repeate the same procedure as before for computing Cos Sim and get the graph
#But first we return to our standard rownames
rownames(mut96_mat_allCDS_norm) <- ssm96_std_rownames
rownames(mut96_mat_allExon_norm) <- ssm96_std_rownames

#For CDS
cos_sim_WGS_WEScds_corr <- cos_sim_matrix(mut96_mat_allWGS_norm, mut96_mat_allCDS_corr)
cos_sim_WGS_WEScds_corr <- as.data.frame(diag(cos_sim_WGS_WEScds_corr))

#For Exon
cos_sim_WGS_WESexon_corr <- cos_sim_matrix(mut96_mat_allWGS_norm, mut96_mat_allExon_corr)
cos_sim_WGS_WESexon_corr <- as.data.frame(diag(cos_sim_WGS_WESexon_corr))

#Plot results
colnames(cos_sim_WGS_WEScds_corr) = "cos_sim"
cos_sim_WGS_WEScds_corr$Tumour_Type = sapply(rownames(cos_sim_WGS_WEScds_corr), function(x) unlist(strsplit(x, "[.]"))[2])

colnames(cos_sim_WGS_WESexon_corr) = "cos_sim"
cos_sim_WGS_WESexon_corr$Tumour_Type = sapply(rownames(cos_sim_WGS_WESexon_corr), function(x) unlist(strsplit(x, "[.]"))[2])

#Let's re-order the cos_sim values based on the median cosine similarity of each tumour type
median_values <- ddply(cos_sim_WGS_WEScds_norm,~Tumour_Type,summarise,median = median(cos_sim,na.rm=TRUE))
order <- median_values$Tumour_Type[order(median_values$median)]
cos_sim_WGS_WEScds_norm$Tumour_Type <- factor(cos_sim_WGS_WEScds_norm$Tumour_Type, levels = order)

cols_pcawg <- pcawg.colour.palette(x = levels(cos_sim_WGS_WEScds_norm$Tumour_Type),scheme = "tumour.subtype")
names(cols_pcawg) <- NULL

boxplot <- ggplot(cos_sim_WGS_WEScds_norm, aes(y=cos_sim, x=Tumour_Type, fill=Tumour_Type)) + 
  geom_boxplot() +
  scale_fill_manual(values=cols_pcawg) +
  ggtitle("Cosine Similarity - SSM96 freq NO-CORRECTED profile \n CDS REGION") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15)))+
  geom_hline(yintercept = 0.9, color="green") +
  geom_hline(yintercept = 0.5, color="red")


pdf("CosSimWGSCDSfreq_NoCorr.pdf",width=12, height=14)
boxplot
dev.off()


# Paired results of Cos Sim values ----------------------------------------

# In this case we wanna plot previous results per tumour type (facet_wrap) and just in CDS or Exon
# The idea is to include lines between paired samples (Corrected or not corrected) to see the effect

# Data we need: cos_sim_WGS_WEScds_norm, cos_sim_WGS_WESexon_norm, 
# cos_sim_WGS_WEScds_corr y cos_sim_WGS_WESexon_corr

# CDS region
cos_sim_WGS_WEScds_corr <- cbind(cos_sim_WGS_WEScds_corr, Corr = rep("YES",dim(cos_sim_WGS_WEScds_corr)[1]),
                                 ID = sapply(rownames(cos_sim_WGS_WEScds_corr), function(x) unlist(strsplit(x,split=".",fixed=TRUE))[1]))
cos_sim_WGS_WEScds_norm <- cbind(cos_sim_WGS_WEScds_norm, Corr = rep("NO",dim(cos_sim_WGS_WEScds_norm)[1]),
                                 ID = sapply(rownames(cos_sim_WGS_WEScds_norm), function(x) unlist(strsplit(x,split=".",fixed=TRUE))[1]))
cos_sim_WGS_WEScds <- rbind(cos_sim_WGS_WEScds_corr,cos_sim_WGS_WEScds_norm)

# Exon region
cos_sim_WGS_WESexon_corr <- cbind(cos_sim_WGS_WESexon_corr, Corr = rep("YES",dim(cos_sim_WGS_WESexon_corr)[1]),
                                 ID = sapply(rownames(cos_sim_WGS_WESexon_corr), function(x) unlist(strsplit(x,split=".",fixed=TRUE))[1]))
cos_sim_WGS_WESexon_norm <- cbind(cos_sim_WGS_WESexon_norm, Corr = rep("NO",dim(cos_sim_WGS_WESexon_norm)[1]),
                                 ID = sapply(rownames(cos_sim_WGS_WESexon_norm), function(x) unlist(strsplit(x,split=".",fixed=TRUE))[1]))
cos_sim_WGS_WESexon <- rbind(cos_sim_WGS_WESexon_corr,cos_sim_WGS_WESexon_norm)

# Plottin!!
boxplot <- ggplot(cos_sim_WGS_WESexon, aes(y=cos_sim, x=Corr, fill=Corr)) + 
  facet_wrap(~ Tumour_Type, scales="free") +
  geom_boxplot() +
  geom_point(aes(x=Corr, y=cos_sim),colour="black",size=2, alpha=0.7) +
  #scale_fill_manual(values=cols_pcawg) +
  ggtitle("Cosine Similarity - SSM96 freq AFTER vs BEFORE k3mer Correction \n EXON region corrected to WGS") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15)))+
  geom_hline(yintercept = 0.9, color="green") +
  geom_hline(yintercept = 0.5, color="red") +
  geom_line(aes(group = ID), size=0.5, colour="grey30", linetype="11", alpha=0.7)


pdf("CosSimWGS_EXON_pairedSamples.pdf",width=18, height=20)
boxplot
dev.off()

# Seleccionamos dos muestras de Uterus-AdenoCA (max y mínima mejora)
Uterus <- cos_sim_WGS_WEScds[cos_sim_WGS_WEScds$Tumour_Type == "Uterus-AdenoCA",]
deltas <- Uterus[Uterus$Corr == "YES","cos_sim"] - Uterus[Uterus$Corr == "NO","cos_sim"]

# Max delta: sample idx 2497 (Delta 0.18061)
# Min delta: sample idx 2482 (Delta -0.19980)

require(MutationalPatterns)
require(SomaticSignatures)

# Percentage of mutations WES in WGS per sample ---------------------------
# HERE WE ARE ALSO INCLUDING THE GENDER DIFFERENCES ANALYSIS - REDUNDANT WITH descriptive...._Gender_Differences.R

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
  facet_wrap( ~ Tumour_Type) +
  geom_boxplot() + scale_y_continuous(trans = "log10") +
  ggtitle("Total number of mutations (SSM96 counts) \n WES CDS data") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12),
                     axis.text.y = element_text(size=14),
                     plot.title = element_text(face = "bold", size = (15))) +
  stat_summary(fun.data = give.n, geom = "text",
               position = position_dodge(width = 0.75),
               size=4, angle=90) + 
  stat_compare_means(label="p.signif") 

pdf("Total_SSMs_CDSwGender_wTest.pdf",width=12, height=14)
boxplot
dev.off()
