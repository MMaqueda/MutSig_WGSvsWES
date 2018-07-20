
###########
# This is for testing the depletion or enrichment of certain SSMs in RoIs
# This is based on the number of counts (SSMs) independently of the trinucleotide seq
###########

# Read the complete mutational profile for WGS, Exon and CDS
WGS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_allWGS.txt",header=TRUE,check.names=FALSE)  
Exon <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/mut96_mat_allExon.txt",header=TRUE,check.names=FALSE) 
CDS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_allCDS.txt",header=TRUE,check.names=FALSE) 

# All three mutational profiles have samples in the same order

# Prepare a data frame where to store all the data. One per RoI (CDS or Exon)

SSMs_CDS_hyper <- as.data.frame(matrix(ncol=7, nrow=dim(WGS)[2]))
colnames(SSMs_CDS_hyper) <- c("Sample_ID", "pval_Depletion", "pval_Enrichment",
                              "adjpval_Depletion", "adjpval_Enrichment",
                              "min_adjpval", "doubling")

SSMs_Exon_hyper <- as.data.frame(matrix(ncol=7, nrow=dim(WGS)[2]))
colnames(SSMs_Exon_hyper) <- c("Sample_ID", "pval_Depletion", "pval_Enrichment",
                               "adjpval_Depletion", "adjpval_Enrichment",
                               "min_adjpval", "doubling")

# Values for the testing

posWGS <- 2861327131
percentageCDS <- 1.23/100
percentageExon <- 4.26/100

# Apply hypergeometric test per sample and store the pvalues

for(i in 1:dim(WGS)[2])
{
  print(paste0("Sample:",i))
  
  # Data for the hypergeometric test with hypothesis in case of depletion (reverse for enrichment):
  # H0: prob of selecting this number of mutations - or fewer - (i.e. CDS region) is no lower than random selection from the WGS
  # H1: prob of selecting this number of mutations - or fewer - (i.e. CDS region) is lower than random selection from the WGS
  
  # For ENRICHMENT, we have to SUBTRACT x BY 1, since if lower.tail is TRUE (default), probabilities are P[X ≤ x], 
  # otherwise, P[X > x]. We subtract x by 1, when P[X ≥ x] is needed. And that's what we wanna test (that number of mutations or higher).
  
  #Depletion test
  pval_CDS_dep <- round(phyper(colSums(CDS)[i], colSums(WGS)[i], posWGS - colSums(WGS)[i], 
                               posWGS * percentageCDS), 4)  
  pval_Exon_dep <- round(phyper(colSums(Exon)[i], colSums(WGS)[i], posWGS - colSums(WGS)[i], 
                                posWGS * percentageExon), 4)  
  
  #Enrichment test
  pval_CDS_enr <- round(phyper(colSums(CDS)[i]-1, colSums(WGS)[i], posWGS - colSums(WGS)[i], posWGS * percentageCDS,
                         lower.tail=FALSE), 4) 
  pval_Exon_enr <- round(phyper(colSums(Exon)[i]-1, colSums(WGS)[i], posWGS - colSums(WGS)[i], posWGS * percentageExon,
                          lower.tail=FALSE), 4)  
  
  # Feed the data frame results: "Sample_ID", "pval_Depletion", "pval_Enrichment", "adjpval_Depletion", 
  # "adjpval_Enrichment", "min_adjpval" and "doubling"
  
  SSMs_CDS_hyper[i,c(1,2,3)] <- c(colnames(WGS)[i],pval_CDS_dep,pval_CDS_enr,
                           round(p.adjust(pval_CDS_dep,"fdr"), 4),
                           round(p.adjust(pval_CDS_enr,"fdr"), 4))
  
  SSMs_Exon_hyper[i,c(1,2,3)] <- c(colnames(WGS)[i],pval_Exon_dep,pval_Exon_enr,
                            round(p.adjust(pval_Exon_dep,"fdr"), 4),
                            round(p.adjust(pval_Exon_enr,"fdr"), 4))
}

# Hacer ajuste por multiple testing en CDS y Exon

SSMs_CDS_hyper[ , "adjpval_Depletion"] <- round(p.adjust(SSMs_CDS_hyper[,"pval_Depletion"],"fdr", 
                                                         n=2*dim(SSMs_CDS_hyper)[1]), 4)
SSMs_CDS_hyper[ , "adjpval_Enrichment"] <- round(p.adjust(SSMs_CDS_hyper[,"pval_Enrichment"],"fdr", 
                                                        n=2*dim(SSMs_CDS_hyper)[1]), 4)

SSMs_Exon_hyper[ , "adjpval_Depletion"] <- round(p.adjust(SSMs_Exon_hyper[,"pval_Depletion"],"fdr", 
                                                          n = 2*dim(SSMs_CDS_hyper)[1]), 4)
SSMs_Exon_hyper[ , "adjpval_Enrichment"] <- round(p.adjust(SSMs_Exon_hyper[,"pval_Enrichment"],"fdr",
                                                           n=2*dim(SSMs_CDS_hyper)[1]), 4)

SSMs_CDS_hyper[ , "min_adjpval"] <- as.numeric(apply(SSMs_CDS_hyper, 1, function(x) 
  min(c(x["adjpval_Depletion"], x["adjpval_Enrichment"]))))

SSMs_CDS_hyper[ , "doubling"] <- 2*SSMs_CDS_hyper[ , "min_adjpval"]
SSMs_CDS_hyper[ , "doubling">1] <- 1

SSMs_Exon_hyper[ , "min_adjpval"] <- min(c(SSMs_Exon_hyper[ , "adjpval_Depletion"], SSMs_Exon_hyper[ , "adjpval_Enrichment"]))
SSMs_Exon_hyper[ , "doubling"] <- 2*SSMs_Exon_hyper[ , "min_adjpval"]


########
save(SSMs_CDS_hyper, SSMs_Exon_hyper,
     file="/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/SSMs_Depletion_Enrichment_RoIs/HypergPERsample.RData")

# Identify samples which show significant difference in CDS/Exon compared to WGS

# In CDS
which(SSMs_CDS_hyper[,"doubling"]<0.05)

# In Exon






if(tumors[i] == "MicroSI")
{
  load("/project/devel/PCAWG/stats_replTime2cluster2drivers.RData")
  MSI_samples <- stats_replTime2cluster2drivers[which(stats_replTime2cluster2drivers[,"msi_riken"] %in% 1),c("sample_id","tumor_type")]
  
  labels <- character()
  for(j in 1:dim(MSI_samples)[1])
  {labels[j] <- paste(MSI_samples$sample_id[j], MSI_samples$tumor_type[j],sep=".")}
  
  CDS <- mut96_mat_CDS[, which(colnames(mut96_mat_CDS) %in% labels)]
  WGS <- mut96_mat_WGS[, which(colnames(mut96_mat_WGS) %in% labels)]
  Exon <- mut96_mat_Exon[,which(colnames(mut96_mat_Exon) %in% labels)]
}
else if(tumors[i] == "Hypermutators") #Parece que falta una muestra de hypermutators (una de colorectal!!)
{
  #De ColoRect van de 2.5M SSMs a 234k SSMs
  #De Uterus es de 280k SSMs
  labels <- c("1582.ColoRect-AdenoCA","2335.ColoRect-AdenoCA","1585.ColoRect-AdenoCA",
              "1572.ColoRect-AdenoCA","1596.ColoRect-AdenoCA","1587.ColoRect-AdenoCA",
              "1615.ColoRect-AdenoCA ","2480.Uterus-AdenoCA")
  
  CDS <- mut96_mat_CDS[, which(colnames(mut96_mat_CDS) %in% labels)]
  WGS <- mut96_mat_WGS[, which(colnames(mut96_mat_WGS) %in% labels)]
  Exon <- mut96_mat_Exon[,which(colnames(mut96_mat_Exon) %in% labels)]
}