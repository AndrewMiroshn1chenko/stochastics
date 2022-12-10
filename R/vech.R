#' Half-Vectorization of a matrix

vech <- function(A, use.Names = TRUE, tol = .Machine$double.eps^0.5) {
  if (is.vector(A)) {
    if (length(A) == 1) {
      return(A)
    } else {
      stop("vech undefined for vectors")
    }
  }

  symm <- is.symmetric.matrix(A, tol)
  if (isFALSE(symm)) {
    stop(sprintf("% must be a numeric and symmetric matrix for half-vectorization."))
  }
  # } else { #if symm is TRUE or NA ...

  full <- vec(A, use.Names)
  stack <- A[lower.tri(A, diag = TRUE)]
  vech.A <- full[stack]

  return(vech.A)
}
