
###########
# This is for testing the depletion or enrichment of certain SSMs in RoIs
# This is based on the number of counts (SSMs) independently of the trinucleotide seq
###########

# Packages to be loaded
require(ggplot2)

# Read the complete mutational profile for WGS, Exon and CDS
mut96_mat_WGS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WGS/txtFiles/mut96_mat_allWGS.txt",header=TRUE,check.names=FALSE)  
mut96_mat_Exon <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/Exon_based/mut96_mat_allExon.txt",header=TRUE,check.names=FALSE) 
mut96_mat_CDS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_allCDS.txt",header=TRUE,check.names=FALSE) 

# All three mutational profiles have samples in the same order

# List of tumour to evaluate

tumors <- c("Breast-DCIS","Breast-AdenoCA","Breast-LobularCA", "Lung-AdenoCA", "ColoRect-AdenoCA", "Prost-AdenoCA","Eso-AdenoCA",
            "CNS-GBM","Head-SCC","Kidney-RCC","Stomach-AdenoCA","Biliary-AdenoCA","Lung-SCC","Skin-Melanoma",
            "Liver-HCC","Thy-AdenoCA","SoftTissue-Liposarc","SoftTissue-Leiomyo","Kidney-ChRCC","CNS-Oligo","Lymph-BNHL",
            "Panc-AdenoCA","Myeloid-AML","Ovary-AdenoCA","CNS-Medullo","CNS-PiloAstro","Uterus-AdenoCA","Bladder-TCC","Cervix-SCC",
            "Panc-Endocrine","Cervix-AdenoCA","Lymph-CLL","Bone-Epith","Bone-Benign","Bone-Osteosarc","Myeloid-MPN","Myeloid-MDS",
            "MicroSI", "Hypermutators")

# Prepare a data frame where to store all the data. One per RoI (CDS or Exon)

SSMs_CDS_hyper <- as.data.frame(matrix(ncol=5,nrow=length(tumors)))
colnames(SSMs_CDS_hyper) <- c("Tumour_Type", "N", "Depleted", "Enriched", "NonSignificant")
SSMs_Exon_hyper <- as.data.frame(matrix(ncol=5,nrow=length(tumors)))
colnames(SSMs_Exon_hyper) <- c("Tumour_Type", "N", "Depleted", "Enriched", "NonSignificant")

# Apply hypergeometric test for each tumor

for(i in 1:length(tumors))
{
  # Filter the specific tumour type
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
  else if(tumors[i] == "Hypermutators")
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
  else #Se trata de un tumor type  
  {
    pattern_tum <- tumors[i] 
    
    WGS <- mut96_mat_WGS[, grep(pattern_tum, colnames(mut96_mat_WGS))]
    Exon <- mut96_mat_Exon[,grep(pattern_tum, colnames(mut96_mat_WGS))]
    CDS <- mut96_mat_CDS[, grep(pattern_tum, colnames(mut96_mat_WGS))]
  }
  
  # Let's merge the data with all SSMs counts (independently of the trinucleotide)
  
  data <- data.frame(SSMs_WGS = colSums(WGS), 
                     SSMs_Exon = colSums(Exon), perc_SSMs_Exon = (colSums(Exon)*100)/colSums(WGS),
                     SSMs_CDS = colSums(CDS), perc_SSMs_CDS = (colSums(CDS)*100)/colSums(WGS))
  data$Sample <- colnames(WGS)
  
  # Prepare data for applying the hypergeometric test
  
  posWGS <- 2861327131
  percentageCDS <- 1.23/100
  percentageExon <- 4.26/100
  numberSamples <- dim(WGS)[2]
  
  # Data for the hypergeometric test with hypothesis in case of depletion (reverse for enrichment):
  # H0: prob of selecting this number of mutations - or fewer - (i.e. CDS region) is no lower than random selection from the WGS
  # H1: prob of selecting this number of mutations - or fewer - (i.e. CDS region) is lower than random selection from the WGS
  
  Data_hyper <- data.frame(x1SSMsCDS = data$SSMs_CDS, x2SSMsExon = data$SSMs_Exon,
                           mSSMsWGS = data$SSMs_WGS,nWGS_woSSMsWGS = posWGS - data$SSMs_WGS,
                           k1CDS = rep(posWGS * percentageCDS,numberSamples), k2Exon = rep(posWGS * percentageExon,numberSamples))
  
  # Apply the hypergeometric test
  
  # For ENRICHMENT, we have to SUBTRACT x BY 1, since if lower.tail is TRUE (default), probabilities are P[X ≤ x], 
  # otherwise, P[X > x]. We subtract x by 1, when P[X ≥ x] is needed. And that's what we wanna test (that number of mutations or higher).
  
  pval_CDS_dep <- apply(Data_hyper, 1, function(x) phyper(x[1],x[3],x[4],x[5]))  #Depletion
  pval_Exon_dep <- apply(Data_hyper, 1, function(x) phyper(x[2],x[3],x[4],x[6])) #Depletion
  
  pval_CDS_enr <- apply(Data_hyper, 1, function(x) phyper(x[1]-1,x[3],x[4],x[5], lower.tail = FALSE)) #Enrichment
  pval_Exon_enr <- apply(Data_hyper, 1, function(x) phyper(x[2]-1,x[3],x[4],x[6],lower.tail = FALSE)) #Enrichment
  
  # Feed the data frame results (("Tumour_Type", "N", "Depleted", "Enriched", "NonSignificant"))
  
  Depleted_CDS <- round((sum(pval_CDS_dep <0.05) / numberSamples)*100, 2)
  Enriched_CDS <- round((sum(pval_CDS_enr <0.05) / numberSamples)*100, 2)
  
  Depleted_Exon <- round((sum(pval_Exon_dep <0.05) / numberSamples)*100, 2)
  Enriched_Exon <- round((sum(pval_Exon_enr <0.05) / numberSamples)*100, 2)
  
  SSMs_CDS_hyper[i, ] <- c(tumors[i], numberSamples,
                           paste0(Depleted_CDS,"%",sep=""),
                           paste0(Enriched_CDS,"%",sep=""),
                           paste0(100- Depleted_CDS-Enriched_CDS,"%",sep=""))
  
  SSMs_Exon_hyper[i, ] <- c(tumors[i], numberSamples,
                           paste0(Depleted_Exon,"%",sep=""),
                           paste0(Enriched_Exon,"%",sep=""),
                           paste0(100- Depleted_Exon- Enriched_Exon,"%",sep=""))
  
  # Generate a plot for that specific tumor type 

  p <- ggplot(data, aes(x=Sample, y=perc_SSMs_Exon)) +
    geom_bar(stat="identity") +
    geom_hline(yintercept = 4.26, color="red",lty="dashed") +
    theme_bw() + theme(axis.text.x = element_blank()) +
    ggtitle(paste0("SSM Perc EXON region in ",tumors[i],"\n Hypergeometric test: ",
                   Depleted_Exon, "% Depleted, ", Enriched_Exon, "% Enriched, Rest % Not Sig",sep="")) 
  
  pdf(paste0("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/SSMs_Depletion_Enrichment_RoIs/Exon_PercSSMs_",tumors[i],".pdf",sep=""))
  plot(p)
  dev.off()
  
  p <- ggplot(data, aes(x=Sample, y=perc_SSMs_CDS)) +
    geom_bar(stat="identity") +
    geom_hline(yintercept = 1.23, color="red",lty="dashed") +
    theme_bw() + theme(axis.text.x = element_blank()) +
    ggtitle(paste0("SSM Perc CDS region in ",tumors[i],"\n Hypergeometric test: ",
                   Depleted_CDS, "% Depleted, ", Enriched_CDS, "% Enriched, Rest % Not Sig",sep="")) 
  
  pdf(paste0("/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/SSMs_Depletion_Enrichment_RoIs/CDS_PercSSMs_",tumors[i],".pdf",sep=""))
  plot(p)
  dev.off()
  
}

save(SSMs_CDS_hyper, SSMs_Exon_hyper,
     file="/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/SSMs_Depletion_Enrichment_RoIs/HypergResults_TOADD.RData")

