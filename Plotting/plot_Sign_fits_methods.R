
#' Compare Contributions from different fits: three methods and seven input data.
#' Plot results per Signature. 

#' @param fits_lsqnonneg  A list with with the fits obtained from fit_to_signatures (SomaticSignatures)
#' @param fits_decSigs  A list with with the fits obtained from deconstructSigs (iterative process)
#' @param fits_QP  A list with with the fits obtained from QP decomposition (SignatureEstimation)
#' @param whichtumour    Character string indicating from which tumour to include in graph title                   
#' @return          Plot
#' REMARK: all the three lists have to have the results (input data) in the same order

plot_Sign_fits_methods <- function(fits_lsqnonneg, fits_decSigs, fits_QP, whichtumour)
{
  # Normalize data in the case of fits_lsqnonneg since it is not norm
  m1 <- lapply(fits_lsqnonneg, function(fit) 
      sweep(as.matrix(fit$contribution), 2, colSums(as.matrix(fit$contribution)), `/`))
  m2 <- fits_decSigs #from deconstructSigs it is already norm
  m3 <- lapply(fits_QP, function(fit) fit$exposure) #We take just the exposure
  
  # Prepare data to plot
  m1_r <- handle_data(m1,"m1")
  m2_r <- handle_data(m2,"m2")
  m3_r <- handle_data(m3,"m3")
  
  # Merge all data
  allresults <- rbind(m1_r,m2_r, m3_r,stringsAsFactors=FALSE)
  
  #And now plot
  
  boxplot <- ggplot(allresults, aes(y=Contribution, x=Data, fill=method),show.legend=F) + 
    facet_wrap(~Signature, scales="free") +
    geom_boxplot() +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, size=10,vjust=0.5),
                       axis.text.y = element_text(size=14),
                       plot.title = element_text(face = "bold", size = (15))) +
    geom_hline(yintercept = 0.05, color="red", lty = "dashed") +
    stat_compare_means(aes(label = ..p.signif.., color="red"),
                       method = "kruskal.test") +
    labs(title=paste("Signature Contribution in",whichtumour, " dataset \n (Kruskal Wallis test between methods)"))
  
  return(boxplot)
}

handle_data <- function(allfits_norm,method)
{
  allresults <- data.frame()
  for(i in 1:length(allfits_norm))
  {
    aux <- as.data.frame(t(allfits_norm[[i]]),stringAsFactors=FALSE)
    aux <- cbind(aux, Data=names(allfits_norm)[i],stringsAsFactors=FALSE)
    allresults <- rbind(allresults, aux, stringsAsFactors=FALSE)
  }
  
  allresults <- melt(allresults)
  colnames(allresults) <- c("Data", "Signature", "Contribution")
  allresults$method <- rep(method, dim(allresults)[1])
  
  return(allresults)
}