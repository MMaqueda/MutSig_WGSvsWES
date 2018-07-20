FrobeniusNorm <- function(M, P, E) {
  sqrt(sum((M-P%*%E)^2))
}

is.wholenumber <- function(x, tol = .Machine$double.eps) {
  abs(x - round(x)) < tol
}
