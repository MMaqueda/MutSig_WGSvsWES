
# Script for generating a data frame with the contribution percentage of all signatures along all tumor types
# Known signatures decomposition (COSMIC). We will only consider results from NNLS (method 1 in the computations)

# We have to obtain the total number of SSMs per tumor type and the number of mutations explained per a
# particular signature in that tumor type to compute the percentage. Contributions less than 0.01 (1%) are not shown.


# List of tumors tested in the work

tumors <- c("Breast-DCIS","Breast-AdenoCA","Breast-LobularCA", "Lung-AdenoCA", "ColoRect-AdenoCA", "Prost-AdenoCA","Eso-AdenoCA",
            "CNS-GBM","Head-SCC","Kidney-RCC","Stomach-AdenoCA","Biliary-AdenoCA","Lung-SCC","Skin-Melanoma",
            "Liver-HCC","Thy-AdenoCA","SoftTissue-Liposarc","SoftTissue-Leiomyo","Kidney-ChRCC","CNS-Oligo","Lymph-BNHL",
            "Panc-AdenoCA","Myeloid-AML","Ovary-AdenoCA","CNS-Medullo","CNS-PiloAstro","Uterus-AdenoCA","Bladder-TCC","Cervix-SCC",
            "Panc-Endocrine","Cervix-AdenoCA","Lymph-CLL","Bone-Epith","Bone-Benign","Bone-Osteosarc","Myeloid-MPN","Myeloid-MDS")

# Create the dataframes

mut_explainedWGS <- as.data.frame(matrix(nrow=30, ncol=length(tumors)))
rownames(mut_explainedWGS) <- paste0("Signature.",seq(1:30))
colnames(mut_explainedWGS) <- tumors

mut_explainedCDS <- as.data.frame(matrix(nrow=30, ncol=length(tumors)))
rownames(mut_explainedCDS) <- paste0("Signature.",seq(1:30))
colnames(mut_explainedCDS) <- tumors

mut_explainedCDScorr <- as.data.frame(matrix(nrow=30, ncol=length(tumors)))
rownames(mut_explainedCDScorr) <- paste0("Signature.",seq(1:30))
colnames(mut_explainedCDScorr) <- tumors


for(i in 1:length(tumors))
{
  # First retrieve the NNLS results for that specific tumor
  # PATH SHOULD BE CHANGED TO WHERE DATA IS STORED
  load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", tumors[i],".RData",sep=""))
  
  # Results from m1 - lsqnonneg - WGS, CDS and CDScorr
  fit_contributionWGS <- fits_m1$WGS$contribution
  fit_contributionCDS <- fits_m1$CDS$contribution
  fit_contributionCDScorr <- fits_m1$CDS_corr$contribution
  
  mut_explainedWGS[,i] <- rowSums(fit_contributionWGS)
  mut_explainedCDS[,i] <- rowSums(fit_contributionCDS)
  mut_explainedCDScorr[,i] <- rowSums(fit_contributionCDScorr)
}

# And now create the percentage data frame from these previous dfs with the mutation counts
# Percentage of mutations explained

mut_expl_perc_CDS <- round((sweep(mut_explainedCDS, 2, colSums(mut_explainedCDS), `/`)*100),1)
mut_expl_perc_WGS <- round((sweep(mut_explainedWGS, 2, colSums(mut_explainedWGS), `/`)*100),1)
mut_expl_perc_CDScorr <- round((sweep(mut_explainedCDScorr, 2, colSums(mut_explainedCDScorr), `/`)*100),1)

mut_expl_perc_CDS <- apply(mut_expl_perc_CDS, c(1,2), function(x) if(x<5) {x <- "-"} else {x<-x})
mut_expl_perc_CDScorr <- apply(mut_expl_perc_CDScorr, c(1,2), function(x) if(x<5) {x <- "-"} else {x<-x})
mut_expl_perc_WGS <- apply(mut_expl_perc_WGS, c(1,2), function(x) if(x<5) {x <- "-"} else {x<-x})

# Optional: save the dataframes
