#' From a mutational matrix (SSM96) with information from different tumour types,
#' only those referring to one or more than one are taken
#' @param Path_mut96_mat  a mutation matrix (SSM96) to be read from a text file. 
#'                        To include path to the specific file
#' @param whichtumour     Tumour type/s to call for the mutation matrix.
#' @param resultsPath     Directory where the results wants to be stored
#' @return                Result is stored in a txt file in indicated resultsPath
#'                        Additionally is returned as a data frame.
#' 
#' #' Example of usage:
#' select_tumour(Path_mut96_mat = "/project/devel/PCAWG/mmaqueda/SSM96Matrix/mut96_mat.txt",
#'     resultsPath = "/project/devel/PCAWG/mmaqueda/SSM96Matrix/",
#'     whichtumour = c("Breast-DCIS","Lung-SCC"))

select_tumour <- function(Path_mut96_mat, resultsPath,whichtumour) {
  
  mut_counts <- read.table(Path_mut96_mat,
                           header=TRUE,check.names=FALSE)
  
  #Independently if there is just one tum type specified or more
  specific_mut_counts <- mut_counts[,unlist(sapply(whichtumour, function(x) grep(x=colnames(mut_counts), pattern= x)))]
  
  write.table(x = specific_mut_counts,
              file = paste0(resultsPath,"mut96_mat_", 
                            paste(whichtumour,collapse="_"),".txt",
                            sep=""),
              sep = "\t")
  
  return(specific_mut_counts)
}