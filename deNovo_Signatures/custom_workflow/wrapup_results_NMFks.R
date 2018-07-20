
#' Function to wrap-up the results from a k's NMF estimation in one list
#' The storage structure is equal to Mut Patterns package for using its plotting functions
#' from the different runs from a NMF estimation such as:
#' estimate <- nmf(mut96_mat, rank=2:14, method="brunet", nrun=100, seed=123456,.options="kvp2")
#' Previous one uses two cores in case they are available. It saves all runs results 
#'
#' @param frobenius    List referring to the output object of kSelectionStableFrob(). 
#'                     This is used to select the best nrun
#' @param nmf_estimation  the NMF.rank object where several k has been estimated and with multiple nrun
#' @return             A list of n elements referring to each k value explored
#'                     Each element is a list of 3: $signatures, $contribution and $reconstructed
#' 
#' Example of usage:
#' nmf_results <- wrapup_results_NMFks(frobenius = estiFrob,
#'     nmf_estimation=estimateLung)
#' 

wrapup_results_NMFks <- function(frobenius,nmf_estimation){
  
  #nrun with best reconstruction error
  best_nrun <- apply(frobenius[[2]],2, which.min)
  
  nmf_results <- list()
  for(i in 1:length(best_nrun))
  {
    nmf_results[[i]] <- list(signatures = basis(nmf_estimation$fit[[i]]@.Data[[best_nrun[i]]]),
                             contribution  = coef(nmf_estimation$fit[[i]]@.Data[[best_nrun[i]]]),
                             reconstructed = basis(nmf_estimation$fit[[i]]@.Data[[best_nrun[i]]]) %*%
                               coef(nmf_estimation$fit[[i]]@.Data[[best_nrun[i]]]))
    
    colnames(nmf_results[[i]]$signatures) <- c(paste("Signature",
                                                     LETTERS[1:(nmf_estimation$measures$rank[i])], sep=""))
    rownames(nmf_results[[i]]$contribution) <- c(paste("Signature",
                                                       LETTERS[1:(nmf_estimation$measures$rank[i])], sep=""))
  }
  return(nmf_results)
}