#' Generating Special Matrices
#' @export I
I <- function(n) {
  # Identity Matrix where number of columns is n.
  D <- diag(1, n)
  dimnames(D) <- list(1:n, 1:n)
  return(D)
}

#' @export J
J <- function(n, m = n) {
  # A matrix of ones with n rows and m columns.
  matrix(c(1), nrow = n, ncol = m, dimnames = list(1:n, 1:m))
}
