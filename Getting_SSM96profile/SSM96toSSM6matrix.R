
#' Transforms/Collapses a SSM96 mutational matrix into the SSM6 base substitutions.
#' @param  Path_mut96_mat a mutation matrix (SSM96) to be read from a text file. 
#'                        Counts! To include path to the specific file.
#' @return                A data frame with the transformed matrix
#' 
#' #' Example of usage:
#' SSM96toSSM6matrix(Path_mut96_mat = "/project/devel/PCAWG/mmaqueda/Lung/mut96_mat_Lung.txt")

SSM96toSSM6matrix <- function(Path_mut96_mat) {
  
  SSM6names <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  mut_counts <- read.table(Path_mut96_mat,
                           header=TRUE,check.names=FALSE)  
  
  SSM6matrix <- rbind(colSums(as.matrix(mut_counts[1:16, ])),
                colSums(as.matrix(mut_counts[17:32, ])),
                colSums(as.matrix(mut_counts[33:48, ])),
                colSums(as.matrix(mut_counts[49:64, ])),
                colSums(as.matrix(mut_counts[65:80, ])),
                colSums(as.matrix(mut_counts[81:96, ])))
  rownames(SSM6matrix) <- SSM6names
  
  return(SSM6matrix)  
}