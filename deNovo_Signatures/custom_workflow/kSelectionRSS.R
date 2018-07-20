
#' Function aiming to choose the best k values based on RSS from a NMF estimation such as:
#' estimate <- nmf(mut96_mat, rank=2:14, method="brunet", nrun=100, seed=123456,.options="kvp2")
#' Previous one uses two cores in case they are available. It saves all runs results 
#'
#' @param nmf_measures    the measures from NMF.rank object where several k has been estimated 
#' @param cutoffRSSrel    The achieved relative percentage of RSS range that is acceptable
#' @return                A vector with one or more k values that could be chosen
#' 
#' Example of usage:
#' kSelectionRSS(nmf_measures = estimate$measures,
#'     cutoffRSSrel=0.85)


kSelectionRSS <- function(nmf_measures,cutoffRSSrel=0.80) {
  
  rangeRSS <- max(nmf_measures$rss) - min(nmf_measures$rss)
  
  #Delta values CCC_i+1 - CCC_i
  diff <- diff(nmf_measures$rss) 
  
  #RSS is always a monotonically decreasing function with k
  selected <- nmf_measures$rank[nmf_measures$rss < (min(nmf_measures$rss) + (rangeRSS * (1-cutoffRSSrel)))]
  return(selected)
}