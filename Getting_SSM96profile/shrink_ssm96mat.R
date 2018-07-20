
#' This function is based on the one developed by Xin Xie (Sept2017)
#' Removes the maximum set of mutation types in a motif matrix that
#' together account for less than or equal to the threshold argument.
#' @param  Path_mut96_mat a mutation matrix (SSM96) to be read from a text file. 
#'                        Counts! To include path to the specific file.
#' @param  threshold      the cutoff cumulative fraction to filter out minor
#'                        types
#' @param resultsPath     Directory where the results wants to be stored
#' @return                Result is stored in a txt file in indicated resultsPath
#'                        Additionally is returned as a data frame.
#' 
#' #' Example of usage:
#' shrink_ssm96mat(Path_mut96_mat = "/project/devel/PCAWG/mmaqueda/SSM96Matrix/mut96_mat.txt",
#'     threshold = 0.01,
#'     resultsPath = "/project/devel/PCAWG/mmaqueda/SSM96Matrix/")

shrink_ssm96mat <- function(Path_mut96_mat, threshold = 0.01, resultsPath) {
  
  mut_counts <- read.table(Path_mut96_mat,
                          header=TRUE,check.names=FALSE)  
  
  #Total number of mutations for all genomes (samples)
  total <- sum(mut_counts)
  
  #Number of counts for each mutation type - Ascend sorted
  sorted_counts_mut_types <- sort(rowSums(mut.counts))
  
  #Detect the minor ones and removed them
  muts_minor <- names(which(cumsum(sorted_counts_mut_types) <= (threshold*total)))
  reduced_mut96_mat <- mut_counts[-(which(rownames(mut_counts) %in% muts_minor)),]
  
  cat("Info: Following mutations types have been removed: \n ", 
      muts_minor)
  
  write.table(x = reduced_mut96_mat,
              file = paste0(resultsPath,"reduced_mut96_mat.txt",
                            sep=""),
              sep = "\t")
  return(reduced_mut96_mat)
}







