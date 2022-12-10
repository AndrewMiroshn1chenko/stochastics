#' Is a matrix symmetric or positive-definite?

# Is matrix square? (A must be a matrix. ).
# Adapted from matrixcalc::is.square.matrix. Argument name changed to be consistent.
is.square.matrix <- function(A) {
  if (!is.numeric(A)) stop(sprintf("% is not a numeric matrix", A))
  if (!is.matrix(A)) stop(sprintf("% is not a matrix", A))
  is.square <- nrow(A) == ncol(A)
  return(is.square)
}

#  @description Is matrix is symmetric? A must be a numeric square matrix with no missing values.
is.symmetric.matrix <- function(A, tol = .Machine$double.eps^0.5) {
  if (anyNA(A)) {
    return(NA)
  }
  if (!is.square.matrix(A)) {
    warning(sprintf("% is not a square matrix", A))
    return(FALSE)
  }

  # Is A and t(A) equal within a tolerance?    # Edited from matrixcalc (PH)
  total.abs <- sum(abs(A - t(A)))
  if (total.abs < tol) { # to avoid being exactly equal
    okay <- TRUE
  } else {
    print("A is not symmetric. Top of the matrix: ")
    print(utils::head(A))
    okay <- FALSE
  }
  # cat("sum( abs(A - t(A)) : ", total.abs, "\n")
  # cat("Total Absolute Difference between Matrix & It's Transpose:", total.abs, "\n")
  return(okay)
}

# Find eigenvalues of a symmetric matrix A with no missing values.
# Paul Hargarten added this function.
find.eval <- function(A, tol = .Machine$double.eps^0.5) {
  # Check if A is symmetric.
  is.symm <- is.symmetric.matrix(A, tol)
  if (is.na(is.symm)) {
    eigenvalues <- NA
  } else if (!is.symm) { # Pass a Negative Number so it is returned false
    # stop(sprintf("%s is not symmetric so no eigenvalues are imputed", A))
    eigenvalues <- -100
  } else if (is.symm) {
    # If A is symmetric, find eigenvalues.
    eigenvalues <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
    # Adjust small eigenvalues to be 0  #(Edited from for loop)
    eigenvalues <- ifelse(abs(eigenvalues) < tol, 0, eigenvalues)
  }
  return(eigenvalues)
}

is.positive.semi.definite <- function(A, tol = .Machine$double.eps^0.5) {
  # Positive semi-definite matrix have non-negative eigenvalues.
  eigenvalues <- find.eval(A, tol)
  if (anyNA(eigenvalues)) {
    return(NA)
  }
  pos.semi <- if (min(eigenvalues) >= 0) TRUE else FALSE
  return(pos.semi)
}

is.positive.definite <- function(A, tol = .Machine$double.eps^0.5) {
  # Positive definite matrix have positive e-values.
  eigenvalues <- find.eval(A, tol)
  pos <- if (anyNA(eigenvalues)) {
    NA
  } else if (min(eigenvalues) > 0) {
    TRUE
  } else {
    FALSE
  }
  return(pos)
}


# Edited check from tmvtnorm so that the symmetric test includes eigenvalues instead of the deteriminant.
# But not used in this function.
checkSymmetricPositiveDefinite <- function(x, name = "sigma") {
  if (!isSymmetric(x, tol = sqrt(.Machine$double.eps))) {
    stop(sprintf("%s must be a symmetric matrix", name))
  }
  if (NROW(x) != NCOL(x)) {
    stop(sprintf("%s must be a square matrix", name))
  }
  if (any(diag(x) <= 0)) {
    stop(sprintf(
      "%s all diagonal elements must be positive",
      name
    ))
  }
  min.eval <- min(eigen(x)$values) # edited from tmvnorm function
  if (det(x) <= 0 & min.eval(x) < 0) {
    stop(sprintf("%s must be positive definite", name))
  }
}
