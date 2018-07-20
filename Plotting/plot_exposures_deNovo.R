
plot_exposures_deNovo <- function (exposures, sampleNames, coord_flip = FALSE, tumor, region) 
{
  # Compute the relative exposures
  contribution <- sweep(exposures, 2, colSums(exposures), `/`)
  
  colnames(contribution) <- sampleNames
  
  SIGNATURES <- paste0("Signature",LETTERS,sep="")
  rownames(contribution) <-  SIGNATURES[1:dim(exposures)[1]]

  m_contribution = melt(contribution)
  colnames(m_contribution) = c("Signature", "Sample", "Contribution")
  
  plot = ggplot(m_contribution, aes(x = factor(Sample), 
                                      y = Contribution, fill = factor(Signature), order = Sample)) + 
    geom_bar(position = "fill", stat = "identity", colour = "black") + 
    scale_fill_discrete(name = "Signature") +
    labs(x = "", y = "Relative contribution") + theme_bw() + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = 6, angle = 90, 
                                     vjust = 0.4)) + 
    theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank()) +
    labs(title = paste0(tumor, " samples ", region, " region"))
  
  if (coord_flip) 
    plot = plot + coord_flip() + xlim(rev(levels(factor(m_contribution$Sample))))
  else plot = plot + xlim(levels(factor(m_contribution$Sample)))
  
  return(plot)
}
