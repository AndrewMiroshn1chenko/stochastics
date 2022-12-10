#' Matrix Trace
tr <- function(A) {
  # if( !is.numeric(A) & !is.matrix(A) )
  #      stop(paste( "A", "must be a numeric matrix."))
  if (nrow(A) != ncol(A)) {
    stop(paste("A", "is not a square matrix"))
  }
  stopifnot(nrow(A) == ncol(A))
  sum(diag(A))
}
