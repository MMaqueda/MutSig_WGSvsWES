
#' Compare Errors from different fits: three methods and seven input data.
#' Plot results per type of input data (each with three methods). 

#' @param fits_lsqnonneg  A list with with the fits obtained from fit_to_signatures (SomaticSignatures)
#' @param fits_decSigs  A list with with the fits obtained from deconstructSigs (iterative process)
#' @param fits_QP  A list with with the fits obtained from QP decomposition (SignatureEstimation)
#' @param mut_mats A list with the corresponding mutational profiles used for previous fitting
#' @param cancer_signatures Corresponding known cancer signatures used for previous fitting
#' @return          A data frame with all the errors. Prepare for plotting
#' REMARK: all objects have to have the results (input data) in the same order

get_errors_diff_methods <- function(fits_lsqnonneg, fits_decSigs, fits_decQP, mut_mats, cancer_signatures)
{
  # Normalize the mutational profiles since they are in counts
  mut_mats <- lapply(mut_mats, function(profile) 
    sweep(as.matrix(profile), 2, colSums(as.matrix(profile)), `/`))
  
  # Get the number of samples for this dataset
  num_samples <- dim(mut_mats[[1]])[2] 
  
  # Errors in method 1 (lsqnonneg)
  # In this case, the exposures should be normalized 
  
  exposures_m1 <- lapply(fits_lsqnonneg, function(fit) 
    sweep(as.matrix(fit$contribution), 2, colSums(as.matrix(fit$contribution)), `/`))
  
  errors_m1 <-lapply(seq(length(fits_lsqnonneg)), function(data) {
    sapply(seq(num_samples), function(i) FrobeniusNorm(mut_mats[[data]][,i],
                                                       as.matrix(cancer_signatures),
                                                       exposures_m1[[data]][,i]))}) 
  names(errors_m1) <- names(fits_lsqnonneg)
  errors_m1 <- as.data.frame(errors_m1)
  errors_m1$SampleID <- colnames(mut_mats[[1]])
  errors_m1 <- melt(errors_m1)
  colnames(errors_m1) <- c("SampleID","Data", "Error")
  errors_m1$method <- rep("m1", dim(errors_m1)[1])

  # Errors in method 2 (deconstructSigs)
  # In this case, the exposures is already normalized
  exposures_m2 <-  fits_decSigs
  errors_m2 <-lapply(seq(length(fits_decSigs)), function(data) {
    sapply(seq(num_samples), function(i) FrobeniusNorm(mut_mats[[data]][,i],
                                                       as.matrix(cancer_signatures),
                                                       exposures_m2[[data]][,i]))}) 
  names(errors_m2) <- names(fits_decSigs)
  errors_m2 <- as.data.frame(errors_m2)
  errors_m2$SampleID <- colnames(mut_mats[[1]])
  errors_m2 <- melt(errors_m2)
  colnames(errors_m2) <- c("SampleID","Data", "Error")
  errors_m2$method <- rep("m2", dim(errors_m2)[1])
  
  # Errors in method 3 (QP)
  # In this case, the exposures are already available
  errors_m3 <-lapply(fits_decQP, function(data) data$errors)
  errors_m3 <- as.data.frame(errors_m3)
  errors_m3$SampleID <- colnames(mut_mats[[1]])
  errors_m3 <- melt(errors_m3)
  colnames(errors_m3) <- c("SampleID","Data", "Error")
  errors_m3$method <- rep("m3", dim(errors_m3)[1])
  
  all_method_errors <- rbind(errors_m1,errors_m2,errors_m3,stringsAsFactors=FALSE)
  
  return(all_method_errors)
}