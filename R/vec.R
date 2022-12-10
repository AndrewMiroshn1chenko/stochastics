#' Stacks a Matrix using matrix operator "vec"

vec <- function(A, use.Names = TRUE) {
  if (is.vector(A)) {
    return(A)
  }
  if (!is.matrix(A)) stop("argument A is not a matrix")
  if (!is.numeric(A)) stop("argument A is not a numeric matrix")
  vec.A <- as.vector(A)

  ## PAUL HARGARTEN ADDED THIS CODE:
  if (use.Names) {
    n <- nrow(A)
    c <- ncol(A)

    # if no dimnames, add # instead.
    if (is.null(dimnames(A)[[1]])) {
      dimnames(A)[[1]] <- 1:n
    }
    if (is.null(dimnames(A)[[2]])) {
      dimnames(A)[[2]] <- 1:c
    }
    names(vec.A) <- paste(rep(dimnames(A)[[1]], c),
      rep(dimnames(A)[[2]], each = n),
      sep = ":"
    )
  }
  return(vec.A)
}
