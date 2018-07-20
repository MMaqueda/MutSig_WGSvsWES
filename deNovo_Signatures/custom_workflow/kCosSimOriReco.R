
#' Function aiming to choose the best k values based on Similarity of Original vs Reconstructed 
#' signal. This comes from a NMF estimation where different k values are explored such as:
#' estimate <- nmf(mut96_mat, rank=2:14, method="brunet", nrun=100, seed=123456,.options="kvp2")
#' Previous one uses two cores in case they are available. It saves all runs results.
#' Similarity is assesses based on Cosinus Similarity metric 
#'
#' @param mut96_mat    the original SSM96 mutational matrix used for computing/estimating NMF
#' @param nmf_results  a list of lists with the wrap-up results for each k 
#'                     This would be the output of function wrapup_results_NMFks()
#' @return             A list of n elements referring to each k value explored
#'                     Each element is a df with samples (rows) and cos sim value (cols)
#' 
#' Example of usage:
#' cosimOriRec_ks <- kCosSimOriReco(mut96_mat = mut96_Lung,
#'     nmf_results= nmf_results)
#' 

kCosSimOriReco <- function(mut96_mat, nmf_results) {
  
  #Cos Sim entre la señal original y reconstruida. 
  #Basado en el package de Mutational Patterns
  
  cos_sim_ori_rec <- lapply(nmf_results, function(x)
  {
    # calculate all pairwise cosine similarities
    cos_sim_ori_rec <- cos_sim_matrix(mut96_mat, x$reconstructed)
    # extract cosine similarities per sample between original and reconstructed
    cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
    
    # Adjust data frame for further plotting with gpplot
    colnames(cos_sim_ori_rec) = "cos_sim"
    cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)
    cos_sim_ori_rec
  }) 
  return(cos_sim_ori_rec)
}