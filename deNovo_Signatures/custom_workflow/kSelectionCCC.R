
#' Function aiming to choose the best k values based on CCC from a NMF estimation such as:
#' estimate <- nmf(mut96_mat, rank=2:14, method="brunet", nrun=100, seed=123456,.options="kvp2")
#' Previous one uses two cores in case they are available. It saves all runs results 
#'
#' @param nmf_measures    the measures from NMF.rank object where several k has been estimated 
#' @param cutoffCCCrel    The relative percentage of CCC range allowed to be decreased
#' @return                A vector with one or more k values that could be chosen
#' 
#' #' Example of usage:
#' kSelectionCCC(nmf_measures = estimate$measures,
#'     cutoffCCCrel=0.85)

kSelectionCCC <- function(nmf_measures,cutoffCCCrel=0.85) {
  
  rangeCCC <- max(nmf_measures$cophenetic) - min(nmf_measures$cophenetic)
  
  #Delta values CCC_i+1 - CCC_i
  diff <- diff(nmf_measures$cophenetic) 
  
  #Rare cases: If CCC is a monotonically increasing or decreasing function
  #Choose those k with CCC above the cutoffCCC(relative) of CCC range.
  if(length(which(diff < 0)) == 0 | length(which(diff > 0)) == 0) {
    selected <- nmf_measures$rank[nmf_measures$cophenetic > (min(nmf_measures$cophenetic) + rangeCCC*cutoffCCCrel)]
  }
  
  #CCC with combination of incremental/decremental groups
  
  else{
    #Detect the decremental points
    decre.points <- which(diff<0)
    for(i in 1:length(decre.points))
    {
      #If next point to drop is below certain margin
      if(nmf_measures$cophenetic[decre.points[i]+1] < (min(nmf_measures$cophenetic) + rangeCCC*cutoffCCCrel)){
        selected <- nmf_measures$rank[1:decre.points[i]]
        break
      } #end if
    } #end for
  } #end else
  
  #Shrink results in case there was/were one or more initial "bad points" - These should not be included
  toremove <- which(nmf_measures$cophenetic[which(nmf_measures$rank %in% selected)] < 
                      (min(nmf_measures$cophenetic) + rangeCCC*cutoffCCCrel))
  
  if(length(toremove) > 0) {
    selected <- selected[-(toremove)]
  }
  return(selected)
}