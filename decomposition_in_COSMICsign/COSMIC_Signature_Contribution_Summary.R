
setwd("/Volumes/cluster/mmaqueda/Results/COSMIC_signatures/Regions_Comparison_AllTumors/")
require(coin)   # This is for using wilcox.exact that handles data with possible ties (mostly our case)
require(reshape)
require(ggplot2)

# Summarize the results of Signature Contributions per region
# A matrix will be constructed per tumor type in order to plot a heatmap showing:
# Value of 1 for sign contribution >0.05 (Wilcox test, pval<5%), 0.5 for >0.01 and 0 for <0.01
# In the case of COSMIC, all of them will have a value of 1.

# List of tumors to analyzed
# Be careful with following tumor types because of the number of samples: Myeloid-MDS, Breast-DCIS and Cervix-AdenoCA

tumors <- c("Breast-DCIS","Breast-AdenoCA","Breast-LobularCA", "Lung-AdenoCA", "ColoRect-AdenoCA", "Prost-AdenoCA","Eso-AdenoCA",
            "CNS-GBM","Head-SCC","Kidney-RCC","Stomach-AdenoCA","Biliary-AdenoCA","Lung-SCC","Skin-Melanoma",
            "Liver-HCC","Thy-AdenoCA","SoftTissue-Liposarc","SoftTissue-Leiomyo","Kidney-ChRCC","CNS-Oligo","Lymph-BNHL",
            "Panc-AdenoCA","Myeloid-AML","Ovary-AdenoCA","CNS-Medullo","CNS-PiloAstro","Uterus-AdenoCA","Bladder-TCC","Cervix-SCC",
            "Panc-Endocrine","Cervix-AdenoCA","Lymph-CLL","Bone-Epith","Bone-Benign","Bone-Osteosarc","Myeloid-MPN","Myeloid-MDS")

#####
# Information directly taken from COSMIC web site
#####

COSMIC <- as.data.frame(matrix(nrow=30,ncol=length(tumors)))
colnames(COSMIC) <- tumors
COSMIC <- apply(COSMIC, c(1,2), function(cell) cell <- 0.0)  # By default, all of them are zeros

# And now retrieve the information from COSMIC
COSMIC[c(1,2,3,5,6,8,10,13,17,18,20,26,30),"Breast-DCIS"] <- 1
COSMIC[c(1,2,3,5,6,8,10,13,17,18,20,26,30),"Breast-AdenoCA"] <- 1
COSMIC[c(1,2,3,5,6,8,10,13,17,18,20,26,30),"Breast-LobularCA"] <- 1
COSMIC[c(1,2,4,5,6,13,17),"Lung-AdenoCA"] <- 1
COSMIC[c(1,5,6,10),"ColoRect-AdenoCA"] <- 1
COSMIC[c(1,5,6),"Prost-AdenoCA"] <- 1
COSMIC[c(1,2,4,5,6,13,17),"Eso-AdenoCA"] <- 1
COSMIC[c(1,5,11),"CNS-GBM"] <- 1
COSMIC[c(1,2,4,5,7,13),"Head-SCC"] <- 1
COSMIC[c(1,5,6,13,27),"Kidney-RCC"] <- 1  # Clear cell (92% RCC) + papillary (8%RCC)
COSMIC[c(1,2,5,13,15,17,18,20,21,26,28),"Stomach-AdenoCA"] <- 1
#COSMIC[c(),"Biliary-AdenoCA"] <- 1  #It's not in the list of COSMIC
COSMIC[c(1,2,4,5,13),"Lung-SCC"] <- 1
COSMIC[c(1,5,7,11,17),"Skin-Melanoma"] <- 1
COSMIC[c(1,4,5,6,12,16,17,22,23,24),"Liver-HCC"] <- 1
COSMIC[c(1,2,5,13),"Thy-AdenoCA"] <- 1
#COSMIC[c(),"SoftTissue-Liposarc"] <- 1  #Is it Chrondosarcoma?
#COSMIC[c(),"SoftTissue-Leiomyo"] <- 1  #Is it Chrondosarcoma?
COSMIC[c(1,5,6),"Kidney-ChRCC"] <- 1  # kidney chromophobe
#COSMIC[c(),"CNS-Oligo"] <- 1  #It's not in the list
#COSMIC[c(),"Lymph-BNHL"] <- 1 #To check
COSMIC[c(1,2,3,5,6,13),"Panc-AdenoCA"] <- 1
COSMIC[c(),"Myeloid-AML"] <- 1  #To check
COSMIC[c(1,3,5),"Ovary-AdenoCA"] <- 1
COSMIC[c(1,5,8),"CNS-Medullo"] <- 1
COSMIC[c(1,5,19),"CNS-PiloAstro"] <- 1
COSMIC[c(1,2,5,6,10,13,14,26),"Uterus-AdenoCA"] <- 1  #An adenocarcinoma is a type of carcinoma
COSMIC[c(1,2,5,10,13),"Bladder-TCC"] <- 1
COSMIC[c(1,2,5,6,10,13,26),"Cervix-SCC"] <- 1
COSMIC[c(1,2,3,5,6,13),"Panc-Endocrine"] <- 1
COSMIC[c(1,2,5,6,10,13,26),"Cervix-AdenoCA"] <- 1
#COSMIC[c(),"Lymph-CLL"] <- 1 #Not sure
#COSMIC[c(),"Bone-Epith"] <- 1
#COSMIC[c(),"Bone-Benign"] <- 1
COSMIC[c(1,2,5,6,13,30),"Bone-Osteosarc"] <- 1
#COSMIC[c(),"Myeloid-MPN"] <- 1  #It's not in the list
#COSMIC[c(),"Myeloid-MDS"] <- 1   #It's not in the list

#####
# Loop for all tumor types
#####

for(i in 1:length(tumors))
{
  # First retrieve the results for that specific tumor
  load(paste0("/Volumes/cluster/mmaqueda/Results/COSMIC_signatures/Methods_Comparison_AllTumors/Methods_COSMIC_", tumors[i],".RData",sep=""))
  
  # Results from m1 - lsqnonneg
  # In this case, contribution must be normalized first!! We can do that by specifying norm=TRUE in 

  #pval01 <- lapply(fits_m1, function(region) compute_pval(region$contribution, 0.01, norm=TRUE))
  #pval05 <- lapply(fits_m1, function(region) compute_pval(region$contribution, 0.05, norm=TRUE))
  
  # Considering results from m2 - deconstructSigs
  
  pval01 <- lapply(fits_m2, function(region) compute_pval(region, 0.01, norm=FALSE))
  pval05 <- lapply(fits_m2, function(region) compute_pval(region, 0.05, norm=FALSE))
  
  # Let's prepare the df to plot
  summary_contr <- as.data.frame(matrix(nrow=30, ncol=6))
  rownames(summary_contr) <- c(paste(rep("Signature.",30), seq(1,30,1), sep=""))
  colnames(summary_contr) <- c("COSMIC","WGS","CDS","CDS_corr","Exon","Exon_corr")
  
  summary_contr <- as.data.frame(apply(summary_contr, c(1,2), function(cell) cell <- 0))  # By default, all of them are zeros
  
  # First column is direct without computations
  summary_contr[,1] <- COSMIC[,i]
  
  # And now change the zero values per 0.5 (>0.01) or 1(>0.05)
  
  for(j in 1:length(pval01))
  {
    for(k in 1:dim(summary_contr)[1])
    {
      if(pval01[[j]][k] < 0.05) #Contribution significantly greater than 0.01
      {summary_contr[k,j+1] <- 0.5}
      
      if(pval05[[j]][k] < 0.05) #Contribution significantly greater than 0.05
      {summary_contr[k,j+1] <- 1}
    }
  }

  # Ready for Plotting
  toplot <- melt(summary_contr)
  toplot$Signature <- rep(c(paste(rep("Signature.",30), seq(1,30,1), sep="")),6)
  
  order <- c(paste(rep("Signature.",30), seq(1,30,1), sep=""))
  toplot$Signature <- factor(toplot$Signature, levels = order)
  
  p <- ggplot(toplot, aes(variable, Signature)) + 
    geom_tile(aes(fill = as.factor(value)), colour="white") +
    scale_fill_manual(values=c("white","grey","black")) +
    ggtitle(paste0("Significant Signature Contribution (Wilcox test) \n",tumors[i])) +
    labs(x = "", y = "") + 
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) + 
    theme(legend.position="none",
          panel.border=element_rect(fill = NA, colour="black",size=0.5),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 330, hjust = 0, colour = "grey50"))

  pdf(paste0("SignatureContribution_m2_",tumors[i],".pdf"))
  plot(p)
  dev.off()
  
}


######
# Aux function
######

compute_pval <- function(contribution,threshold, norm)
{
  # Normalize contribution if it is in counts (m1 is in counts not percentage)
  
  if(norm==TRUE)
    contribution <- sweep(as.matrix(contribution), 2, colSums(as.matrix(contribution)), `/`)
  
  pval <- apply(contribution, 1,  function(signature) {
    test <- wilcox.exact(signature, mu=threshold, alternative = "greater")
    return(test$p.value)})
  
  return(pval)
}
