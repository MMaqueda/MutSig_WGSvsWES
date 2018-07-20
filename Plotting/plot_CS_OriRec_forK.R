
#' This function plots a barplot over all samples showing the Cos Sim
#' between the original and the reconstructed signal for that specific k
#' 
#' @param CSOriReco_ks A list with as elements as k values estimated 
#'                     each element is a df with cos_sim value
#' @return  Directly the plot 

#######

#For instance the input would be the output of....:
#skin_CSOriReco_ks <- kCosSimOriReco(skin_melanoma, skin_NMFks_results)

plot_CS_OriRec_forK <- function(CSOriReco_ks,kvalue,nmf_measures){
  
  kcase <- CSOriReco_ks[[which(nmf_measures$rank == kvalue)]]
  
  plot <- ggplot(kcase, aes(y=cos_sim, x=sample)) +
    geom_bar(stat="identity", fill = "skyblue4") +
    coord_cartesian(ylim=c(0.8, 1)) +
    # coord_flip(ylim=c(0.8,1)) +
    ylab("Cosine similarity\n original VS reconstructed") +
    xlab("") +
    theme_bw() +
    theme(panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(), strip.background = element_blank(),
          axis.text.x=element_text(size=5,angle=90)) +
    # Add cut.off line for cos sim
    geom_hline(aes(yintercept=.95))
  
  return(plot)
}
