
Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
source("/project/devel/PCAWG/mmaqueda/Rscripts/pcawg.colour.palette.R")
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggplot2))

#' This is an adaptation from MutationalPatterns package (plot spectrum)
#' But adapted to SSM96 profile for whole set of samples: mean and sd
#' If samples want to be plotted separately use 'plot_profile_SSM96_normalized.R'
#' Plot the relative or absolute counts (argument relative)
#'
#' @param mutSSM96 The standard mutational SSM96 matrix
#' @param wSD Boolean to indicate whether to include error bars in the plot or not
#' @return  A plot object (ggplot2) ready to plot

plot_spectrum_SSM96 <- function (mutSSM96, legend = TRUE, relative=F, wSD=TRUE, ymax=0.2, metric = "mean", tumor, region) 
{
  COLORS6 = c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE")
  
  SUBSTITUTIONS = c('C>A','C>G','C>T','T>A','T>C','T>G')
  
  C_TRIPLETS = c(
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT")
  
  T_TRIPLETS = c(
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT")
  
  CONTEXTS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))
  
  colors = COLORS6
  context = CONTEXTS_96
  substitution = rep(SUBSTITUTIONS, each = 16)
  substring(context, 2, 2) = "."
  
  df = data.frame(substitution = substitution, context = context)
  
  # Relative contribution if specified
  if(relative==T) 
    mutSSM96 <- as.data.frame(mutSSM96/colSums(mutSSM96))

  rownames(mutSSM96) = NULL
  
  # Mean and sd
  mean <-apply(mutSSM96, 1, function(SSM96type)  mean(SSM96type, na.rm=TRUE))
  sd <-apply(mutSSM96, 1, function(SSM96type)  sd(SSM96type, na.rm=TRUE))
  median <-apply(mutSSM96, 1, function(SSM96type)  median(SSM96type, na.rm=TRUE))
  
  # Gather all information for plotting
  if (metric == "mean")
  {df2 = cbind(df, mean)}
  else # metric = median
  {df2 = cbind(df, median)}
    
  df3 = melt(df2, id.vars = c("substitution", "context"))
  
  # Plotting
    plot = ggplot(data = df3, aes(x = context, y = value, 
                                  fill = substitution, width = 1)) + 
    geom_bar(stat = "identity", colour = "black", size = 0.2) + scale_fill_manual(values = colors) + 
    facet_grid(variable ~ substitution) + ylab("Relative contribution") + 
    guides(fill = FALSE) + theme_bw() + 
    theme(axis.title.y = element_text(size = 12, vjust = 1), 
            axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
            axis.text.x = element_text(size = 5, angle = 90, 
                                       vjust = 0.4), strip.text.x = element_text(size = 9), 
            strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
            panel.spacing.x = unit(0, "lines")) +
    labs(title = paste0(tumor, " samples", region, " region"))
    
    if(wSD == TRUE)
      {plot <- plot + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.05)}
    
    if(!is.na(ymax))
    {plot <- plot + coord_cartesian(ylim = c(0, ymax)) + 
      scale_y_continuous(breaks = seq(0, ymax, 0.1))}
    
  
  return(plot)
  
}
