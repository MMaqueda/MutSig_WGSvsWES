#' Normalizes a mutational matrix initially expressed in counts per mutation type
#' and returns it into frequences per mutation type
#' @param Path_mut96_mat      a mutation matrix (SSM96) to be read from a text file. 
#'                        To include path to the specific file
#' @param resultsPath     Directory where the results wants to be stored
#' @return                Result is stored in a txt file in indicated resultsPath
#'                        Additionally is returned as a data frame.
#' 
#' #' Example of usage:
#' normalize_ssm96mat(Path_mut96_mat = "/project/devel/PCAWG/mmaqueda/SSM96Matrix/mut96_mat.txt",
#'     resultsPath = "/project/devel/PCAWG/mmaqueda/SSM96Matrix/")

normalize_ssm96mat <- function(Path_mut96_mat, resultsPath) {
  
  mut_counts <- read.table(Path_mut96_mat,
                           header=TRUE,check.names=FALSE)
  
  frequences_mut96_mat <- sweep(mut_counts, 2, colSums(mut_counts), `/`)
  
  write.table(x = frequences_mut96_mat,
              file = paste0(resultsPath,"freqnorm_mut96_mat.txt",
                            sep=""),
              sep = "\t")
  return(frequences_mut96_mat)
}