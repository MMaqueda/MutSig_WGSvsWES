
Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
source("/project/devel/PCAWG/mmaqueda/Rscripts/pcawg.colour.palette.R")
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggplot2))

# Constants
kNtBase <- c("A", "C", "G", "T")
kMutBase <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
kNumOfMutCtx3 <- 96

#' This function plots the heatmap of all samples based on their SSM96 profile
#' It is based on Xin Xie corresponding functions
#' First, an auxiliary function to prepare the data
#' @param mut96_mat The mutational SSM96 matrix
#' @return  A plot object (ggplot2) ready to plot


# Transforms a motif matrix for heatmap plot on 96 mutation types using ggplot2.

manipulateXheatmap <- function(mut96mat, log10 = TRUE) {
  
  data <- sweep(mut96mat, 2, colSums(mut96mat), `/`)
  
  checkNAs <- apply(data, 2, function(x) all(is.na(x)))
  data[,which(checkNAs == TRUE)] <- 0
  
  data <- apply(data,c(1,2),function(x)
  {
    if(x<0.0001) log10(0.01)
    else log10(x*100) 
  })
  
  data <- as.data.frame(data)
  
  sample.names <- colnames(data)
  n.samples <- length(sample.names)
  
  data <- stack(data)
  data <- as.data.frame(data$values)
  colnames(data) <- "value"
  data$mut_base <- rep(rep(kMutBase, each = 4 * 4), times = n.samples)
  data$fwd_base <- rep(rep(kNtBase, each = 4), times = 6 * n.samples)
  data$nxt_base <- rep(kNtBase, times = 6 * 4 * n.samples)
  data$sample_name <- rep(sample.names, each = kNumOfMutCtx3)
  return(data)
}


plot_heatmap_SSM96 <- function(mut96mat) {
  
  data <- manipulateXheatmap(mut96mat)
  
  plot <- ggplot(data, aes(x = nxt_base, y = fwd_base, fill = value)) +
    facet_grid(sample_name ~ mut_base, labeller = label_context) +
    geom_tile() +
    scale_y_discrete(limits = rev(kNtBase)) +
    scale_fill_gradient2(name = "Percentage (log10)", limits= c(-2, 2),
                         low = "white", mid = "gold", high = "red",
                         midpoint = 0,
                         guide = guide_colourbar(title.position = "top",
                                                 barwidth = 4)) +
    #coord_equal() + #For squares - doable for 20-30 samples maximum
    coord_fixed(ratio = 0.3) +
    theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"),
          strip.text.y=element_text(size=2,angle = 0),
          strip.text.x=element_text(size=4,angle = 0),
          legend.position="bottom",
          axis.text.y=element_text(size=1),
          axis.text.x=element_text(size=2),
          axis.ticks = element_blank()) +
    labs(x = "Three Prime Base", y = "Five Prime Base") 
  
  return(plot)
}