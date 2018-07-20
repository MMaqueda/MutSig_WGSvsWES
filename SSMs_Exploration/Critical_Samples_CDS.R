
# Table with the number of samples with critical number of SSMs. 
# Compute the number of samples <100, <50, <25 and all zero per tumor type.
# ONLY APPLIED TO CDS! (at least a priori)

mut96_mat_CDS <- read.table("/project/devel/PCAWG/mmaqueda/Data/SSM96Matrix_WES/CDS_based/mut96_mat_allCDS.txt",header=TRUE,check.names=FALSE) 

# List of tumour to evaluate

tumors <- c("Breast-DCIS","Breast-AdenoCA","Breast-LobularCA", "Lung-AdenoCA", "ColoRect-AdenoCA", "Prost-AdenoCA","Eso-AdenoCA",
            "CNS-GBM","Head-SCC","Kidney-RCC","Stomach-AdenoCA","Biliary-AdenoCA","Lung-SCC","Skin-Melanoma",
            "Liver-HCC","Thy-AdenoCA","SoftTissue-Liposarc","SoftTissue-Leiomyo","Kidney-ChRCC","CNS-Oligo","Lymph-BNHL",
            "Panc-AdenoCA","Myeloid-AML","Ovary-AdenoCA","CNS-Medullo","CNS-PiloAstro","Uterus-AdenoCA","Bladder-TCC","Cervix-SCC",
            "Panc-Endocrine","Cervix-AdenoCA","Lymph-CLL","Bone-Epith","Bone-Benign","Bone-Osteosarc","Myeloid-MPN","Myeloid-MDS",
            "MicroSI", "Hypermutators")

# Prepare a data frame where to store all the data. 

SSMs_CDS_critical <- as.data.frame(matrix(ncol=6,nrow=length(tumors)))
colnames(SSMs_CDS_critical) <- c("Tumour_Type", "N", "<100", "<50", "<25", "NonSSMs")

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
  }
  else if(tumors[i] == "Hypermutators")
  {
    #De ColoRect van de 2.5M SSMs a 234k SSMs
    #De Uterus es de 280k SSMs

   labels <- c("1582.ColoRect-AdenoCA","2335.ColoRect-AdenoCA","1585.ColoRect-AdenoCA",
                "1572.ColoRect-AdenoCA","1596.ColoRect-AdenoCA","1587.ColoRect-AdenoCA",
                "1615.ColoRect-AdenoCA ","2480.Uterus-AdenoCA")
    CDS <- mut96_mat_CDS[, which(colnames(mut96_mat_CDS) %in% labels)]
  }
  else #Se trata de un tumor type  
  {
    pattern_tum <- tumors[i] 
    CDS <- mut96_mat_CDS[, grep(pattern_tum, colnames(mut96_mat_CDS))]
  }

  
  # Let's sum the SSMs per sample
  totSSMs <- colSums(CDS)

  # Prepare data for applying the hypergeometric test
  numberSamples <- dim(CDS)[2]
  
  # Feed the data frame results 
  SSMs_CDS_critical[i, ] <- c(tumors[i], numberSamples, length(which(totSSMs < 100)),
                           length(which(totSSMs < 50)),length(which(totSSMs < 25)),
                           length(which(totSSMs == 0)))
}

save(SSMs_CDS_critical,
     file="/project/devel/PCAWG/mmaqueda/Data/Comparison_SSM96_WGS_WES/SSMs_CDS_critical.RData")
