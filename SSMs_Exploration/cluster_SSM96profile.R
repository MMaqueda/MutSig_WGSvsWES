
#' Computes the hierarchical clustering of samples based on their SSM96 profile
#' The distance metric is the Cosine Similarity between samples
#' This function is based on the cluster_signatures() from MutPatt package
#' It can be used for samples or signatures

#' @param profiles_matrix     a mutation matrix SSM96 (samples or signatures) 
#' @param method     Method to apply for the hierarchical clustering
#' @return                a hclust object


cluster_SSM96profile <- function (profiles_matrix, method = "complete") 
{
  profiles_matrix <- profiles_matrix[,colSums(profiles_matrix) > 0]
  
  sim = cos_sim_matrix(profiles_matrix, profiles_matrix)  #cos_sim_matrix is from MutationalPattern
  dist = as.dist(1 - sim)
  hc_cos = hclust(dist, method = method)
  return(hc_cos)
}