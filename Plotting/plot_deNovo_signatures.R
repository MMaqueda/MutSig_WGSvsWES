
Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
source("/project/devel/PCAWG/mmaqueda/Rscripts/pcawg.colour.palette.R")
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggplot2))

#' This is to plot the Signatures obtained from deNovo computation
#' They will be plot with SD values (from iterations)
#' 
#' @param processes This is the matrix Processes obtained from deNovo computations
#'                  which includes the signature contribution
#' @param processesSD This is the matrix ProcessesStd obtained from deNovo comp
#'                    which includes the SD from the different iterations
#' @param tumour  Character with the name of the tumor type related (for title purposes)
#' @param region  Character with the name of the region related (for title purposes)                     
#' @return  A plot object (ggplot2) ready to plot

plot_deNovo_signatures <- function(processes, processesSD, legend = TRUE, ymax=0.2, tumor, region) 
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
  
  SIGNATURES <- paste0("Signature",LETTERS,sep="")
  
  # Prepare signatures contribution
  contribution <- as.data.frame(processes)
  colnames(contribution) <- SIGNATURES[1:dim(contribution)[2]]
  rownames(contribution) <- NULL
  
  # Gather all information for plotting
  df2 = cbind(df, contribution)
  df3 = melt(df2, id.vars = c("substitution", "context"))
  colnames(df3) <- c("substitution", "Context","Signature", "Contribution")
  
    # Add the SDs values for signature contribution
  df3$SD <- melt(as.data.frame(processesSD))$value
  
  # Plotting
  plot = ggplot(data = df3, aes(x = Context, y = Contribution, 
                                fill = substitution, width = 1)) + 
    geom_bar(stat = "identity", colour = "black", size = 0.2) + scale_fill_manual(values = colors) + 
    facet_grid(Signature ~ substitution) + ylab("Relative contribution") + 
    coord_cartesian(ylim = c(0, ymax)) + 
    scale_y_continuous(breaks = seq(0, ymax, 0.1)) + 
    geom_errorbar(aes(ymin = Contribution - SD, ymax = Contribution + SD), width = 0.05) +
    guides(fill = FALSE) + theme_bw() + 
    theme(axis.title.y = element_text(size = 12, vjust = 1), 
          axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
          axis.text.x = element_text(size = 5, angle = 90, 
                                     vjust = 0.4), strip.text.x = element_text(size = 9), 
          strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
          panel.spacing.x = unit(0, "lines")) +
    labs(title = paste0(tumor, " samples", region, " region"))
  
  return(plot)
  
}
