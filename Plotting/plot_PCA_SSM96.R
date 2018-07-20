
Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
source("/project/devel/PCAWG/mmaqueda/Rscripts/pcawg.colour.palette.R")
source("/project/devel/PCAWG/mmaqueda/Rscripts/normalize_ssm96mat.R")
suppressPackageStartupMessages(library(pls))

#' This function plots the PCA over all samples based on their SSM96 profile
#' @param mut96_mat The mutational SSM96 matrix: FREQ not COUNTS
#' @return  Directly plots (scoreplot) the PCA. First 4PC by default

plot_PCA_SSM96 <- function(mut96_mat){
  
  pca <-  prcomp(t(mut96_mat),scale=TRUE)
  
  #summary(pca) #If we wanna assess the prop of var % per PC...
  
  Tumour_Type <- sapply(colnames(mut96_mat), function(x) unlist(strsplit(x, "[.]"))[2])
  idx <- sapply(colnames(mut96_mat), function(x) unlist(strsplit(x, "[.]"))[1])
  
  cols_pcawg <- pcawg.colour.palette(x = Tumour_Type,
                                     scheme = "tumour.subtype")
  
  #The colors assignated to Lung SCC and AdenoCA are white or light grey so they are not seen..
  # color.map <- function(Tumour_Type) {if (Tumour_Type=="Lung-SCC") "chocolate2" 
  # else if (Tumour_Type=="Lung-AdenoCA") "aquamarine3" }
  # runnerpoint <- unlist(lapply(Tumour_Type, color.map))
  
  return(scoreplot(pca,pch=16, main="PCA", col=cols_pcawg,
                   comps=1:4,labels=idx,pretty.xlabels))
}





