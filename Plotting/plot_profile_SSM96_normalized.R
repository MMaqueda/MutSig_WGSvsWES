
plot_profile_SSM96_normalized <- function (norm_mut_matrix, ymax = 0.2, condensed = FALSE) 
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
  rownames(norm_mut_matrix) = NULL
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  value = NULL
  if (condensed) {
    plot = ggplot(data = df3, aes(x = context, y = value, 
                                  fill = substitution, width = 1)) + geom_bar(stat = "identity", 
                                                                              colour = "black", size = 0.2) + scale_fill_manual(values = colors) + 
      facet_grid(variable ~ substitution) + ylab("Relative contribution") + 
      coord_cartesian(ylim = c(0, ymax)) + scale_y_continuous(breaks = seq(0, 
                                                                           ymax, 0.1)) + guides(fill = FALSE) + theme_bw() + 
      theme(axis.title.y = element_text(size = 12, vjust = 1), 
            axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
            axis.text.x = element_text(size = 5, angle = 90, 
                                       vjust = 0.4), strip.text.x = element_text(size = 9), 
            strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
            panel.spacing.x = unit(0, "lines"))
  }
  else {
    plot = ggplot(data = df3, aes(x = context, y = value, 
                                  fill = substitution, width = 0.6)) + geom_bar(stat = "identity", 
                                                                                colour = "black", size = 0.2) + scale_fill_manual(values = colors) + 
      facet_grid(variable ~ substitution) + ylab("Relative contribution") + 
      coord_cartesian(ylim = c(0, ymax)) + scale_y_continuous(breaks = seq(0, 
                                                                           ymax, 0.1)) + guides(fill = FALSE) + theme_bw() + 
      theme(axis.title.y = element_text(size = 12, vjust = 1), 
            axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
            axis.text.x = element_text(size = 5, angle = 90, 
                                       vjust = 0.4), strip.text.x = element_text(size = 9), 
            strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank())
  }
  return(plot)
}