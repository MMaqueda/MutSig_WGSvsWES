
#' Function aiming to choose the best k values based on Similarity of Original vs Reconstructed 
#' signal. This comes from a NMF estimation where different k values are explored such as:
#' estimate <- nmf(mut96_mat, rank=2:14, method="brunet", nrun=100, seed=123456,.options="kvp2")
#' Previous one uses two cores in case they are available. It saves all runs results.
#' Similarity is assessed based on Cosinus Similarity metric. This summmarizes the results from
#' kCosSimOriReco() by computing the mean and sd per k value
#'
#' @param CosSim    a list of dataframes with the CosSim values per sample per k
#'                  this object can be the output of kCosSimOriReco() function
#' @param nmf_estimation  The estimation object from nmf()
#' @return             A dataframe with the mean and sd of cos sim per rank value explored
#' 
#' Example of usage:
#' summCosSim <- summaryCosSim(CosSim = cosimOriRec_ks,
#'     nmf_estimation= estimateLung)
 
summaryCosSim <- function(CosSim,nmf_estimation){
  
  mean <- unlist(lapply(CosSim, function(x) mean(x$cos_sim)))
  sd <- unlist(lapply(CosSim, function(x) sd(x$cos_sim)))
  
  summary <- data.frame(mean=mean, sd=sd, 
                        rank=paste("Rank",nmf_estimation$measures$rank,sep=""))
  return(summary)
}
