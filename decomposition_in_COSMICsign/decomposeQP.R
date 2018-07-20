#' decomposeQP Function
#'
#' This function allows to get the optimal solution by using dual method to solve the quadratic programming problem.
#' @param m observed turmor profile vector for a single patient/sample, 96 by 1. m is normalized.
#' @param P signature profile matrix, 96 by N(N = # signatures, COSMIC: N=30)
#' @param control some control parameter that can be passed into the solve.QP function
#' @keywords quadratic programming
#' @export
#' @examples
#' decomposeQP(tumorBRCA[,1], signaturesCOSMIC)

decomposeQP <- function(m, P, ...){
  # N: how many signatures are selected
  N = ncol(P)
  # G: matrix appearing in the quatric programming objective function
  G = t(P) %*% P
  # C: matrix constraints under which we want to minimize the quatric programming objective function.
  C <- cbind(rep(1,N), diag(N))
  # b: vector containing the values of b_0.
  b <- c(1,rep(0,N))
  # d: vector appearing in the quatric programming objective function
  d <- t(m) %*% P

  #Solve quadratic programming problem
  out = quadprog::solve.QP(Dmat = G, dvec = d, Amat = C, bvec = b, meq = 1)

  #Some exposure values are negative, but very close to 0
  #Change these neagtive values to zero and renormalized
  exposures = out$solution
  exposures[exposures < 0] = 0
  exposures = exposures/sum(exposures)

  # return the exposures
  return(exposures)
}

