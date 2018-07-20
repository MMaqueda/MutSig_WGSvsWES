
#' Function aiming to choose the best k values based on the stability of Frobenius norm values
#' from the different runs from a NMF estimation such as:
#' estimate <- nmf(mut96_mat, rank=2:14, method="brunet", nrun=100, seed=123456,.options="kvp2")
#' Previous one uses two cores in case they are available. It saves all runs results 
#'
#' @param mut96_mat    the original mutational SSM96 matrix usef for computing/estimating NMF
#' @param nmf_estimation  the NMF.rank object where several k has been estimated and with multiple nrun
#' @return             A list with two elements: the first element with a vector with one or more k values
#'                     the second with a df with frob norm per nrun (rows) per k value (cols)
#' 
#' Example of usage:
#' estiFrob <- kSelectionStableFrob(mut96_mat = mut96_Lung,
#'     nmf_estimation=estimateLung)
#' 
#' And to inspect results.....         
#' estiFrob[[1]] #selected
#' boxplot(estiFrob[[2]],use.cols=TRUE)  #plot Frob norm values


kSelectionStableFrob <- function(mut96_mat, nmf_estimation) {
  
  reconstruct <- lapply(nmf_estimation$fit, function(k) {lapply(k@.Data, function(nrun) basis(nrun) %*% coef(nrun))})
  differences <- lapply(reconstruct, function(k) {lapply(k, function(nrun) mut96_mat - nrun)})
  frobenius <- lapply(differences, function(k) {unlist(lapply(k, function(nrun) norm(as.matrix(nrun), type="f")))})
  frobenius <- as.data.frame(frobenius)
  colnames(frobenius) <- paste("rank",nmf_estimation$measures$rank,sep="")
  
  #Select those ones with a range 10 times below than its mean
  
  means <- apply(frobenius, 2, mean)
  range <- apply(frobenius, 2, function(x) max(x)-min(x))
  selected <- nmf_estimation$measures$rank[which(range < means/10)]
  
  return(list(selected,frobenius))
}