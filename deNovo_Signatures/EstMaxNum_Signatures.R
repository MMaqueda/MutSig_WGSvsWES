
#' This function is adapted from the one developed by Xin Xie (September 2017).
#' Estimates the maximum number of possible signatures based on the given number
#' of genomes. When deciphering signatures, a range of possible numbers need to
#' be tested. This range normally starts from 2 and to such estimated maximum
#' number that follows an exponential-fit. The fitted model is built upon six
#' paired points taken from Alexandrov L.B et al Cell Reports (2013).
#' @param  n.genomes      the number of genomes for signature construction
#' @return                the estimated maximum number of possible signatures

EstMaxNumOfSignatures <- function(n.genomes) {
  
  #These values are taken from Alexandrov L.B et al (FigS1B,Fig3C)
  kNumOfGenomes <- c(10, 20, 30, 50, 100, 200)
  kMaxNumOfSignatures <- c(3, 5, 7, 10, 15, 20)
  
  fit.data <- data.frame(n_genome = kNumOfGenomes, n_sig = kMaxNumOfSignatures)
  model <- nls(n_sig ~ a + b * log(n_genome), data = fit.data,
               start = list(a = 0, b = 0))
  estimate <- ceiling(predict(model, data.frame(n_genome = n.genomes)))
  
  return(estimate)
}

