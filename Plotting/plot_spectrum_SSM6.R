
Sys.setenv(BIOCINSTALLER_ONLINE_DCF = F)
source("/project/devel/PCAWG/mmaqueda/Rscripts/pcawg.colour.palette.R")
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(ggplot2))

#' This is an adaptation from MutationalPatterns package (plot spectrum)
#' Plot spectrum SSM6 per sample, per tumour type or all (argument by)
#' Plot the relative or absolute counts (argument relative)
#'
#' @param mutSSM6 The mutational SSM6 matrix
#' @return  A plot object (ggplot2) ready to plot


plot_spectrum_SSM6 <- function (mutSSM6, plotby=c("ALL","TUMOUR","SAMPLE"), 
                                legend = TRUE, relative=F) 
{
  #Fixed colors to be used for each base sustitution
  colors = c("#2EBAED", "#000000", "#DE1C14",
              "#D4D2D2", "#ADCC54", "#F0D0CE")
  
  # #Information in the form of samples x mutations
  tmutSSM6 <- t(mutSSM6) 
  # 
  # #Relative contribution if specified
  if(relative==T) 
    df2 <- as.data.frame(tmutSSM6/rowSums(tmutSSM6))
  else df2 <- as.data.frame(tmutSSM6) 
  # 
  #Grouping variable to use
  if (plotby == "ALL") {by = "all"}
  else if (plotby == "TUMOUR") {by <- sapply(rownames(tmutSSM6), function(x) unlist(strsplit(x, "[.]"))[2]) }
  else {by=rownames(tmutSSM6)} #Per sample
  # 
  # Add by info to df
  df2$by <- by 
   
  # Reshape
  df3 <- melt(df2, id.vars = "by")
   
  # Count number of mutations per mutation type
  counts <- melt(tmutSSM6, measure.vars = colnames(tmutSSM6))
  df4 <- cbind(df3, counts$value)
  colnames(df4)[4] = "nmuts"
   
  # Calculate the mean and the stdev on the value for each group broken
  # down by type + variable
   x <-ddply(df4, c("by", "variable"), summarise, mean = mean(value,na.rm = TRUE), 
             stdev = sd(value,na.rm = TRUE))

  info_x <- ddply(df4, c("by"), summarise, total_individuals = sum(value),
                 total_mutations = sum(nmuts))

  x <- merge(x, info_x)

  info_type <- data.frame(sub_type = c("C>A", "C>G", "C>T",
                                      "T>A", "T>C", "T>G"),
                         variable = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))

  x <- merge(x, info_type)

  x$total_mutations = prettyNum(x$total_mutations, big.mark = ",")
  x$total_mutations = paste("No. mutations =", as.character(x$total_mutations))

  # Define positioning of error bars
  x$error_pos = x$mean
  
  #Make barplot
  plot = ggplot(data = x, aes(x = sub_type, y = mean, fill = variable, group = sub_type)) +
    geom_bar(stat = "identity") +
    #coord_flip() +
    scale_fill_manual(values = colors, name = "Point mutation type") +
    theme_bw() + xlab("") +
    scale_y_continuous(limits = c(0, 1)) +
    ylab("Relative Contribution") +
    theme(axis.text.x = element_blank(),panel.grid.major.x = element_blank(),
          axis.ticks = element_blank())

  # check if standard deviation error bars can be plotted
  if (sum(is.na(x$stdev)) > 0)  {warning("No standard dev error bars can be plotted: there is only one sample per mutation spectrum")}
  else plot = plot + geom_errorbar(aes(ymin = error_pos - stdev, ymax = error_pos + stdev), width = 0.2)
  
  # Facetting
  if (plotby=="ALL") {plot = plot + facet_wrap(~total_mutations)}
  else if (plotby == "TUMOUR") {plot = plot + facet_wrap(by ~ total_mutations)}
  else {
    # my_grob = grobTree(textGrob(paste("Test",seq(1:9)), 
    #                             x=0.8,  y=0.95, hjust=0,
    #                             gp=gpar(col="blue", fontsize=7, fontface="italic")))
    plot = plot + facet_wrap( ~ by, scales="free") + 
    theme(strip.background = element_blank(), strip.text=element_text(size=3),
          panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0, "lines"),
          axis.text.x=element_text(size=7,angle=90,vjust=0.5,hjust=1),
          panel.grid.minor = element_blank(),panel.grid.major = element_blank()) 
    #annotation_custom(my_grob)
    }
  
  # Legend
  if (legend == FALSE)
    plot = plot + theme(legend.position = "none")
  return(plot)
  }
