
#' Compare Contributions from different fits and plot them per Signature. Assuming that input data 
#' has been obtained from fit_to_signatures function

#' @param list_of_fits    A list with with the fits obtained from fit_to_signatures. Each list element
#'                        corresponds to a fit with different data origin 
#' @param whichtumour    Character string indicating from which tumour to include in graph title                   
#' @return          Plot
#' 


plot_multiple_Sign_fits <- function(allfits, whichtumour, norm, decQP=FALSE)
{
  # Normalize data since it is in counts (from fit_to_signatures)
  if(norm==TRUE){
    allfits_norm <- lapply(allfits, function(fit) 
      sweep(as.matrix(fit$contribution), 2, colSums(as.matrix(fit$contribution)), `/`))
  }
  else allfits_norm <- allfits #from deconstructSigs it is already norm
  
  if(decQP==TRUE) {allfits_norm  <- lapply(allfits, function(fit) fit$exposure)}
  
  # Prepare data to plot
  
  allresults <- data.frame()
  for(i in 1:length(allfits_norm))
  {
    aux <- as.data.frame(t(allfits_norm[[i]]),stringAsFactors=FALSE)
    aux <- cbind(aux, Data=names(allfits_norm)[i],stringsAsFactors=FALSE)
    allresults <- rbind(allresults, aux, stringsAsFactors=FALSE)
  }
  
  # Plot the data + Kruskal Wallis test (compare all groups)
  
  allresults <- melt(allresults)
  colnames(allresults) <- c("Data", "Signature", "Contribution")
  
  allresults$Data <- factor(allresults$Data, levels = names(allfits))
  
  boxplot <- ggplot(allresults, aes(y=Contribution, x=Data, fill=Data),show.legend=F) + 
    facet_wrap(~Signature, scales="free") +
    geom_boxplot() +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, size=10,vjust=0.5),
                       axis.text.y = element_text(size=14),
                       plot.title = element_text(face = "bold", size = (15))) +
    geom_hline(yintercept = 0.05, color="red", lty = "dashed") +
    stat_compare_means(aes(label = ..p.signif..),
                       method = "wilcox.test", ref.group = "WGS", paired=TRUE) +
    labs(title=paste("Signature Contribution in",whichtumour, " dataset (Paired Wilcoxon test - WGS ref)"))

  return(boxplot)
}