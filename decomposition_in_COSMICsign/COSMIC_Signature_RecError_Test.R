
# This script is for applying a t.test for testing if there are significant differences between CDS and CDS corr
# or Exon and Exon corr in terms of reconstruction error.
# We will only consider the results from m1 method -NNLS- (since there were no significant differences between 
# the three methods, NNLS, heuristic and QP. A data frame is generated with the results.
# Additionally, we will generate a data frame where the mean +/- sd reconstruction error is stored

# List of tumors to analyzed
# Be careful with following tumor types because of the number of samples: Myeloid-MDS, Breast-DCIS and Cervix-AdenoCA

tumors <- c("Breast-DCIS","Breast-AdenoCA","Breast-LobularCA", "Lung-AdenoCA", "ColoRect-AdenoCA", "Prost-AdenoCA","Eso-AdenoCA",
            "CNS-GBM","Head-SCC","Kidney-RCC","Stomach-AdenoCA","Biliary-AdenoCA","Lung-SCC","Skin-Melanoma",
            "Liver-HCC","Thy-AdenoCA","SoftTissue-Liposarc","SoftTissue-Leiomyo","Kidney-ChRCC","CNS-Oligo","Lymph-BNHL",
            "Panc-AdenoCA","Myeloid-AML","Ovary-AdenoCA","CNS-Medullo","CNS-PiloAstro","Uterus-AdenoCA","Bladder-TCC","Cervix-SCC",
            "Panc-Endocrine","Cervix-AdenoCA","Lymph-CLL","Bone-Epith","Bone-Benign","Bone-Osteosarc","Myeloid-MPN","Myeloid-MDS")

#Create data frames

summary_tests <- as.data.frame(matrix(nrow=length(tumors), ncol=2))
rownames(summary_tests) <- tumors
colnames(summary_tests) <- c("CDSvsCDScorr","ExonvsExoncorr")

# We'll take profit and generate a data frame where the mean +/- sd reconstruction error is stored

recerror_avg <- as.data.frame(matrix(nrow=length(tumors), ncol=5))
rownames(recerror_avg) <- tumors
colnames(recerror_avg) <- c("WGS","CDS","corrected CDS", "Exon", "corrected Exon")

means <- as.data.frame(matrix(nrow=length(tumors), ncol=5))
rownames(means) <- tumors
colnames(means) <- c("WGS","CDS","corrected CDS", "Exon", "corrected Exon")

#####
# Loop for all tumor types
#####

for(i in 1:length(tumors))
{
  # First retrieve the results for that specific tumor
  # CHANGE THE PATH TO WHERE RESULTS ARE STORES
  load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", tumors[i],".RData",sep=""))
  
  # Results from m1 - lsqnonneg
  errors <- errors_methods[errors_methods$method == "m1" , ] 
  
  recerror_avg[i,] <- paste0(round(tapply(errors$Error, errors$Data, function(x) mean(x, na.rm=TRUE)),3), " +/- ", 
                             round(tapply(errors$Error, errors$Data, function(x) sd(x, na.rm=TRUE)),3), sep="")
  
  means[i,] <- tapply(errors$Error, errors$Data, function(x) mean(x, na.rm=TRUE))
  
  # CDS vs CDS_corr

   pvalCDS <- t.test(x = errors[errors$Data == "CDS", "Error"],
                     y = errors[errors$Data == "CDS_corr", "Error"],
                     alternative = "two.sided",
                     paired = TRUE)
  
   pvalExon <- t.test(x = errors[errors$Data == "Exon", "Error"],
                      y = errors[errors$Data == "Exon_corr", "Error"],
                      alternative = "two.sided",
                      paired = TRUE)
  
   # Store the results
   summary_tests[i,1] <- round(pvalCDS$p.value,5)
  summary_tests[i,2] <- round(pvalExon$p.value,5)
  
}

summary_tests$CDSvsCDScorr <- p.adjust(summary_tests$CDSvsCDScorr, method="fdr")
summary_tests$ExonvsExoncorr <- p.adjust(summary_tests$ExonvsExoncorr, method="fdr")

# Let's fix (prepare them for publishing) the p-value results for the report

fix_adjpval <- function(pval)
{
  if(pval<0.001) {
    pval <- "<.001"}else if(pval<0.05) {
      pval <- "<.05"}else if(pval<0.01) {
        pval <- "<.01"}else{pval <- "NS"}
  
  return(pval)
}

summary_tests <- apply(summary_tests, c(1,2), function(x) fix_adjpval(x))
recerror_avg <- cbind(recerror_avg[,c(1,2,3)], summary_tests[,1],recerror_avg[,c(4,5)], summary_tests[,2])

#Optional: save the dataframes with the results
