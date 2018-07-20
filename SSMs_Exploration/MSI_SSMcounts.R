
# MSI Samples. Number of SSMs in RoI vs WGS -------------------------------

# Read the mutational profiles for the different regions
mut96_mat_WGS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_allWGS.txt",header=TRUE,check.names=FALSE)  
mut96_mat_Exon <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/mut96_mat_allExon.txt",header=TRUE,check.names=FALSE) 
mut96_mat_CDS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_allCDS.txt",header=TRUE,check.names=FALSE) 

# Select the MSI samples

load("/project/devel/PCAWG/stats_replTime2cluster2drivers.RData")
MSI_samples <- stats_replTime2cluster2drivers[which(stats_replTime2cluster2drivers[,"msi_riken"] %in% 1),c("sample_id","tumor_type")]

labels <- character()
for(i in 1:dim(MSI_samples)[1])
{labels[i] <- paste(MSI_samples$sample_id[i], MSI_samples$tumor_type[i],sep=".")}

MSI_WGS <- mut96_mat_WGS[, which(colnames(mut96_mat_WGS) %in% labels)]
MSI_Exon <- mut96_mat_Exon[, which(colnames(mut96_mat_Exon) %in% labels)]
MSI_CDS <- mut96_mat_CDS[, which(colnames(mut96_mat_CDS) %in% labels)]

# Let's merge the data with all SSMs counts (independently of the trinucleotide)

MSI_SSMs <- data.frame(SSMs_WGS = colSums(MSI_WGS), 
                       SSMs_Exon = colSums(MSI_Exon), perc_SSMs_Exon = (colSums(MSI_Exon)*100)/colSums(MSI_WGS),
                       SSMs_CDS = colSums(MSI_CDS), perc_SSMs_CDS = (colSums(MSI_CDS)*100)/colSums(MSI_WGS))
rownames(MSI_SSMs) <- colnames(MSI_WGS)

# Apply the hypergeometric test on those for DEPLETION

data <- MSI_SSMs

posWGS <- 2861327131
percentageCDS <- 1.23/100
percentageExon <- 4.26/100
#numberSamples <- dim(MSI_SSMs)[1]
numberSamples <- dim(Eso_SSMs)[1]
numberSamples <- dim(Kidney_SSMs)[1]

Data_hyper_kid <- data.frame(x1SSMsCDS = data$SSMs_CDS -1,
                         x2SSMsExon = data$SSMs_Exon -1,
                         mSSMsWGS = data$SSMs_WGS,
                         nWGS_woSSMsWGS = posWGS - data$SSMs_WGS,
                         k1CDS = rep(posWGS * percentageCDS,numberSamples), 
                         k2Exon = rep(posWGS * percentageExon,numberSamples))

# Hypothesis
# H0: prob of selecting this number of mutations (i.e. CDS region) is no lower than random selection from the WGS
# H1: prob of selecting this number of mutations (i.e. CDS region) is lower than random selection from the WGS

pval_kid_CDS <- apply(Data_hyper_kid, 1, function(x) phyper(x[1],x[3],x[4],x[5]))  #Depletion
pval_kid_Exon <- apply(Data_hyper_kid, 1, function(x) phyper(x[2],x[3],x[4],x[6])) #Depletion

pval_kid_CDS_enrich <- apply(Data_hyper_kid, 1, function(x) phyper(x[1],x[3],x[4],x[5],
                                                                       lower.tail = FALSE)) #Enrichment
pval_kid_Exon_enrich <- apply(Data_hyper_kid, 1, function(x) phyper(x[2],x[3],x[4],x[6],
                                                                       lower.tail = FALSE)) #Enrichment

# Let's put all this data in the same place

MSI_SSMs <- cbind(MSI_SSMs, pval_Exon = pval_Exon, pval_CDS = pval_CDS,
                  pval_Exon_enrich = pval_Exon_enrich, pval_CDS_enrich = pval_CDS_enrich)

# Prepare data for plotting
MSI_SSMs$Sample <- rownames(MSI_SSMs)

pval_Dep <- formatC(MSI_SSMs$pval_CDS,format = "e", digits = 2)
pval_Enr <- formatC(MSI_SSMs$pval_CDS_enrich,format = "e", digits = 2)

labels <- sapply(seq(1,16), function(p) paste("Depleted ",pval_Dep[p]," Enriched ",pval_Enr[p],
                                              "  SSMs_WGS=",MSI_SSMs$SSMs_WGS[p],sep=""))

# Plotting

#sum(MSI_SSMs$pval_CDS_enrich <0.05)
  
pdf("MSI_CDS.pdf")
ggplot(MSI_SSMs, aes(x=Sample, y=perc_SSMs_CDS)) +
  geom_bar(stat="identity",fill="lightgrey") +
  geom_text(aes(label=labels), angle=90, color="black", size=2.2 , hjust=1) +
  geom_hline(yintercept = 1.23, color="red",lty="dashed") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=10)) +
  ggtitle("MSI samples: SSM Percentage (Hypergeometric Depleted or Enriched P-value) \n CDS region (No multiple correction)") +
  annotate("label", x = 2, y = 2, label = "Depleted 5/16 \n Enriched 9/16")
dev.off()


# Eshopagus samples -------------------------------------------------------

# It's supposed that in Esophagus case, mutations tend to be outside the exon

# Select the Eso samples

pattern_Eso <- "Eso-AdenoCA" 

Eso_WGS <- mut96_mat_WGS[, grep(pattern_Eso, colnames(mut96_mat_WGS))]
Eso_Exon <- mut96_mat_Exon[,grep(pattern_Eso, colnames(mut96_mat_WGS))]
Eso_CDS <- mut96_mat_CDS[, grep(pattern_Eso, colnames(mut96_mat_WGS))]

# Let's merge the data with all SSMs counts (independently of the trinucleotide)

Eso_SSMs <- data.frame(SSMs_WGS = colSums(Eso_WGS), 
                       SSMs_Exon = colSums(Eso_Exon), perc_SSMs_Exon = (colSums(Eso_Exon)*100)/colSums(Eso_WGS),
                       SSMs_CDS = colSums(Eso_CDS), perc_SSMs_CDS = (colSums(Eso_CDS)*100)/colSums(Eso_WGS))
rownames(Eso_SSMs) <- colnames(Eso_WGS)
data <- Eso_SSMs

# Hyper already computed

Eso_SSMs <- cbind(Eso_SSMs, pval_Exon = pval_eso_Exon, pval_CDS = pval_eso_CDS,
                  pval_Exon_enrich = pval_eso_Exon_enrich, pval_CDS_enrich = pval_eso_CDS_enrich)

# Prepare data for plotting
Eso_SSMs$Sample <- rownames(Eso_SSMs)

pval_Dep <- formatC(Eso_SSMs$pval_Exon, format = "e", digits = 2)
pval_Enr <- formatC(Eso_SSMs$pval_Exon_enrich,format = "e", digits = 2)

labels <- sapply(seq(1,length(pval_Dep)), function(p) paste("Depleted ",pval_Dep[p]," Enriched ",pval_Enr[p],
                                              "  SSMs_WGS=",Eso_SSMs$SSMs_WGS[p],sep=""))

# Plotting

#sum(Eso_SSMs$pval_Exon<0.05)

pdf("Eso_Exon.pdf")
ggplot(Eso_SSMs, aes(x=Sample, y=perc_SSMs_Exon)) +
  geom_bar(stat="identity",fill="lightgrey") +
  geom_text(aes(label=labels), angle=90, color="black", size=1, hjust=0.5) +
  geom_hline(yintercept = 4.26, color="red",lty="dashed") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=5)) +
  ggtitle("Esophagus samples: SSM Percentage (Hypergeometric Depleted or Enriched) \n EXON region (No multiple correction)") +
  annotate("label", x = 10, y = 7.5, label = "Depleted 94/97 \n Enriched 1/97")
dev.off()

# Theoretical minimum for Enrichment --------------------------------------

# How much of the CDS should be mutated to be enriched? Just exploring the MSI samples. Data in Data_hyper
# Instead of phyper, we should use qhyper.


qval_CDS_enrich <- apply(Data_hyper, 1, function(x) qhyper(0.05,x[3],x[4],x[5],
                                                                   lower.tail = FALSE)) #Enrichment

toplot <- data.frame(sample = rownames(MSI_SSMs), cutoff_SSM_Enrich = qval_CDS_enrich, actual_SSMs = MSI_SSMs$SSMs_CDS)
toplot <- melt(toplot)

pdf("MSI_Cutoff_Enrichment.pdf")
ggplot(toplot, aes(y=value, x=sample, fill= variable)) + 
  geom_bar(stat="identity",position=position_dodge()) +
  ggtitle("ENRICHMENT Hypergeometric test - SSMs counts \n CDS REGION") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, size=12)) 
dev.off()


