#' findSigExposures Function
#' wrapper function
#' This function allows to obtain the optimal solution by specifying quadratic programming or simulated annealing to solve the optimization problem.
#' @param M observed turmor profile matrix for all the patient/sample, 96 by G. G is the number of patients. Each column can be mutation
#' counts, or mutation probabilities. Each column will be normalized to sum up to 1.
#' @param P signature profile matrix, 96 by N(N = # signatures, COSMIC: N=30)
#' @param decompostion.method which method is selected to get the optimal solution: decomposeQP or decomposeSA
#' @keywords optimal method: QP or SA
#' @export
#' @examples
#' E1 = findSigExposures(tumorBRCA, signaturesCOSMIC, decomposeQP)
#' sigsBRCA = c(1,2,3,5,6,8,13,17,18,20,26,30)
#' E2 = findSigExposures(tumorBRCA, signaturesCOSMIC[, sigsBRCA], decomposeQP)
#' E3 = findSigExposures(tumorBRCA[, 1:10], signaturesCOSMIC, decomposeSA, list(maxit=1000, temperature=100))
#' E4 = findSigExposures(tumorBRCA[, 1:10], signaturesCOSMIC[, sigsBRCA], decomposeSA, list(maxit=2000))
#' E5 = findSigExposures(round(tumorBRCA*10000), signaturesCOSMIC, decomposeQP)

findSigExposures <- function(M, P, decomposition.method = decomposeQP, norm=TRUE, ...) {
  ## process and check function parameters
  ## M, P
  M = as.matrix(M)
  P = as.matrix(P)
  if(nrow(M) != nrow(P))
    stop("Matrices 'M' and 'P' must have the same number of rows (mutations types).")
  if(any(rownames(M) != rownames(P)))
    stop("Matrices 'M' and 'P' must have the same row names (mutations types).")
  if(ncol(P) == 1)
    stop("Matrices 'P' must have at least 2 columns (signatures).")
  ## decomposition.method
  if(!is.function(decomposition.method))
    stop("Parameter 'decomposition.method' must be a function.")
  ## normalize M by column (just in case it is not normalized)
  if(norm==TRUE){M = apply(M, 2, function(x){x/sum(x)})}
  
  
  ## find solutions
  ## matrix of signature exposures per sample/patient (column)
  exposures = apply(M, 2, decomposition.method, P, ...)
  rownames(exposures) = colnames(P)
  colnames(exposures) = colnames(M)
  
  ## compute estimation error for each sample/patient (Frobenius norm)
  errors = sapply(seq(ncol(M)), function(i) FrobeniusNorm(M[,i],P,exposures[,i]))
  names(errors) = colnames(M)
  
  return(list(exposures=exposures, errors=errors))
}