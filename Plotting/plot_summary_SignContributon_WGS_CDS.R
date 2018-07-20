
#' Boxplot of the signature contribution from WGS, CDS and CDScorr for a specific tumor type
#' Kruskal test to find significant differences between them

#' @param SignContributionWGS    a df with the contributions for WGS
#' @param SignContributionCDS    a df with the contributions for CDES
#' @param SignContributionCDScorr    a df with the contributions for CDS corrected
#' @param whichtumour    Character string indicating from which tumour to include in graph title
#' @param norm    Boolean. Indicate if dataframe must be normalized (TRUE)  - if contribution is 
#'                in counts instead of fraction (i.e. output fit_to_signatures -NNLS) or not (FALSE) 
#'                (i.e. output of whichSignatures - heuristic)  
#' @return          Plot
#' 

# Example:
# Contributions have to be retrieved e.g. from method 1 (NNLS)
#load(paste0("/project/devel/PCAWG/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", "Lung-AdenoCA",".RData",sep=""))
#plot_summary_SignaturesContribution_WGS_CDS(fits_m1$WGS$contribution, fits_m1$CDS$contribution,fits_m1$CDS_corr$contribution,"Lung-AdenoCA",norm=TRUE)

plot_summary_SignaturesContribution_WGS_CDS <- function(SignContributionWGS, SignContributionCDS, SignContributionCDScorr, whichtumour, norm){
  
  # Normalize data if it is in counts instead of percentage (this is the case for method 1 - NNLS)
  if(norm==TRUE)
  {
    SignContributionWGS <- sweep(as.matrix(SignContributionWGS), 2, colSums(as.matrix(SignContributionWGS)), `/`)
    SignContributionCDS <- sweep(as.matrix(SignContributionCDS), 2, colSums(as.matrix(SignContributionCDS)), `/`)
    SignContributionCDScorr <- sweep(as.matrix(SignContributionCDScorr), 2, colSums(as.matrix(SignContributionCDScorr)), `/`)
    }
  else
    {
      SignContributionWGS <- SignContributionWGS
      SignContributionCDS <- SignContributionCDS
      SignContributionCDScorr <- SignContributionCDScorr
    }
  
  # Set up all the colors possible
  all.sigs <- c(paste(rep("Signature.",30), seq(1,30,1), sep=""))              
  # all.colors  <- c("#023FA5","#BEC1D4","#D6BCC0","#BB7784","gold","#4A6FE3","#8595E1","#B5BBE3",
  #                  "#E6AFB9","#E07B91","#D33F6A","#11C638","#8DD593","#C6DEC7","#EAD3C6","#F0B98D",
  #                  "#EF9708","#0FCFC0","#9CDED6","#D5EAE7","#F3E1EB","#F6C4E1","#F79CD4","#866097",
  #                  "#008941","#A30059","#008080","#8B0000","#F4A460","#663399")
  
  df1 <- melt(SignContributionWGS)
  colnames(df1) <- c("Signature","SampleID","Contribution")
  df1$Region <- rep("WGS", dim(df1)[2])
  
  df2 <- melt(SignContributionCDS)
  colnames(df2) <- c("Signature","SampleID","Contribution")
  df2$Region <- rep("CDS", dim(df2)[2])
  
  df3 <- melt(SignContributionCDScorr)
  colnames(df3) <- c("Signature","SampleID","Contribution")
  df3$Region <- rep("CDScorr", dim(df3)[2])
  
  df4 <- rbind(df1,df2,df3)
  
  #To keep the order
  df4$Signature <- factor(df4$Signature, levels = all.sigs) # to keep order
  
  boxplot <- ggplot(df4, aes(y=Contribution, x=Signature, fill=Region)) + 
    geom_boxplot() +
    ylim(0,1) +
    ggtitle(paste("Signature Contribution in",whichtumour, " dataset")) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, size=10,vjust=0.5),
                       axis.text.y = element_text(size=14),
                       plot.title = element_text(face = "bold", size = (15)))+
    stat_compare_means(aes(group = Region), label = "p.signif", hide.ns=TRUE,method="kruskal.test") + 
    geom_hline(yintercept = 0.05, color="red", lty = "dashed") 
  
  return(boxplot)
}