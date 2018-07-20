
#' Plot the mean and sd values for Cos Similarity between Original and Reconstructed
#' samples per all k values estimated in NMF  
#'
#' @param CosSim_summary    a dataframe with the mean and sd values per k assessed
#'                          Typically the output from summaryCosSim() function
#' @return          Plot
#' 

#Input: the list of dataframes from previous function 


plot_summaryCosSim_OrRec <- function(CosSim_summary){
  
  plot <- ggplot(CosSim_summary, aes(y=mean, x=rank)) +
      geom_bar(stat="identity", fill = "chocolate2") +
      coord_cartesian(ylim=c(min(CosSim_summary$mean - CosSim_summary$sd), 
                             max(CosSim_summary$mean + CosSim_summary$sd))) +
      ylab("Avg Cosine similarity\n original VS reconstructed") +
      xlab("") +
      theme_bw() +
      theme(panel.grid.minor.y=element_blank(),
            panel.grid.major.y=element_blank(), strip.background = element_blank(),
            axis.text.x=element_text(size=10,angle=90)) +
      #Adding the std dev bar
      geom_errorbar(aes(ymin = CosSim_summary$mean - CosSim_summary$sd, 
                        ymax = CosSim_summary$mean + CosSim_summary$sd), width = 0.2)+
      # Add cut.off line for cos sim
      geom_hline(aes(yintercept=.98),col="red")
    
  return(plot)
}
