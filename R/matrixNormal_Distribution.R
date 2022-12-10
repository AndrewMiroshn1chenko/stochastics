#' The Matrix Normal Distribution
dmatnorm <- function(
  A,
                     M,
                     U,
                     V,
                     tol = .Machine$double.eps^0.5,
                     log = TRUE
                     ) {
  n <- nrow(A)
  p <- ncol(A)

  # Checks
  if (is.data.frame(A)) A <- as.matrix(A)
  if (sum(dim(A) == dim(M)) != 2) stop("M must have same dimensions as A.")
  check_matnorm(s = 1, M, U, V, tol)

  # The Log Density
  log.dens <- (-n * p / 2) * log(2 * pi) - p / 2 * log(det(U)) - n / 2 * log(det(V)) +
    -1 / 2 * tr(solve(U) %*% (A - M) %*% solve(V) %*% t(A - M))

  # Return
  if (log) {
    return(log.dens)
  } else {
    return(exp(log.dens))
  }
}

# #'@importFrom LaplacesDemon logdet() .  This is the Function used:
logdet <- function(x) {
  2 * sum(log(diag(chol(x))))
}
# Uses logdet(U) instead of log(det(U)) which could calculate .
# The log(det(U)) and log(det(V)) terms have same terms but are positive.

# Uses logdet function. Any difference using LaplaceNormal?
dmatnorm.logdet <- function(A, M, U, V,
                            tol = .Machine$double.eps^0.5,
                            log = TRUE) {
  n <- nrow(A)
  p <- ncol(A)

  # Checks
  if (is.data.frame(A)) A <- as.matrix(A)
  if (sum(dim(A) == dim(M)) != 2) stop("M must have same dimensions as A.")
  check_matnorm(s = 1, M, U, V, tol)

  # The Log Density
  log.dens <- (-n * p / 2) * log(2 * pi) - p / 2 * logdet(U)
  -n / 2 * logdet(V) +
    -1 / 2 * tr(solve(U) %*% (A - M) %*% solve(V) %*% t(A - M))

  # Return
  if (log) {
    return(log.dens)
  } else {
    return(exp(log.dens))
  }
}



#' @rdname matrixNormal_Distribution
#' @param Lower	 The n x p matrix of lower limits for CDF.
#' @param Upper	 The n x p matrix of upper limits for CDF.
#' @inheritParams mvtnorm::pmvnorm
pmatnorm <- function(
  Lower = -Inf,
                     Upper = Inf,
                     M,
                     U,
                     V,
                     tol = .Machine$double.eps^0.5,
                     keepAttr = TRUE,
                     algorithm = mvtnorm::GenzBretz(),
                     ...
                     ) {
  if (utils::packageVersion("mvtnorm") < "1.1-2") {
    warning("New argument added to `mvtnorm v. 1.1-2`. Please upgrade to avoid error when passing `keepAttr`.")
  }

  n <- nrow(M)
  p <- ncol(M)

  # Checks
  check_matnorm(s = 1, M, U, V, tol)

  # Convert the matrices to lower
  if (is.matrix(Lower)) {
    lower <- vec(Lower)
  } else {
    if (is.vector(Lower) & Lower == -Inf) {
      lower <- -Inf
    } else {
      stop("The lower limit must be a numeric matrix or -Inf.")
    }
  }

  if (is.matrix(Upper)) {
    upper <- vec(Upper)
  } else {
    if (is.vector(Upper) & Upper == Inf) {
      upper <- Inf
    } else {
      stop("The upper limit must be a numeric matrix or Inf.")
    }
  }
  # Calculating the probablity
  prob <- mvtnorm::pmvnorm(
    lower, upper,
    mean = vec(M),
    corr = NULL,
    sigma = kronecker(U, V),
    algorithm = algorithm,
    ...,
    keepAttr = keepAttr
  )
  warning("The covariance matrix is standardized. ")

  return(prob)
}

#' @param s The number of observations desired to simulate from the matrix normal. Defaults to 1. Currently has no effect but acts as a placeholder in future releases.
# #'@inheritParams mvtnorm::rmvnorm
#' @param method String specifying the matrix decomposition used to determine the matrix root of the Kronecker product of U and V in \code{rmatnorm}. Possible methods are eigenvalue decomposition ("eigen"), singular value decomposition ("svd"), and Cholesky decomposition ("chol"). The Cholesky (the default) is typically fastest, but not by much though. Passed to \code{\link[mvtnorm]{rmvnorm}}.

#' @import mvtnorm
#' @export rmatnorm

rmatnorm <- function(
  s = 1,
                     M,
                     U,
                     V,
                     tol = .Machine$double.eps^0.5,
                     method = "chol"
  ) {
  if (utils::packageVersion("mvtnorm") < "1.1-2") {
    warning("New argument added to `mvtnorm v. 1.1-2`. Please upgrade to avoid error.")
  }

  # Convert all to matrices -- added 5/2/20
  M <- as.matrix(M)
  U <- as.matrix(U)
  V <- as.matrix(V)

  n <- nrow(M)
  p <- ncol(M)

  # Checks
  # (Symmetry is already checked)
  check_matnorm(s, M, U, V, tol)

  # Vectorizing and sampling from rmvnorm
  if (utils::packageVersion("matrixNormal") <= "0.0.5") {
    warning("The construction of sigma has been found to be incorrect. Please upgrade to new version.")
  }
  # Sigma <- kronecker(U,V) #incorrect -- thanks @prockenschaub, https://github.com/phargarten2/matrixNormal/issues/1
  Sigma <- kronecker(V, U)

  vec.X <- mvtnorm::rmvnorm(1, vec(M), Sigma, method = method, checkSymmetry = FALSE)
  # pre0.9_9994 = FALSE #uses later version of package.

  # Reputting back into a matrix
  X <- matrix(vec.X,
    nrow = n, ncol = p,
    dimnames = list(rownames(U), colnames(V))
  )

  return(X)
}

# cov(vec(A))  #should be 1

# Check to make sure the parameters in MatrixNormal match.
check_matnorm <- function(s, M, U, V, tol) {
  if (!(s > 0)) stop("s must be > 0. s = ", s, call. = FALSE)
  if (anyNA(M)) {
    stop("M contains missing values.", call. = FALSE)
  }
  if (anyNA(U)) {
    stop("U contains missing values.")
  }
  if (anyNA(V)) {
    stop("V contains missing values.")
  }
  if (nrow(M) != nrow(U)) {
    stop("The mean matrix M has different sample size than the scale sample size
         matrix U. M has ", dim(M)[[1]], "rows, and U has ", dim(U)[[1]], ".")
  }
  if (ncol(M) != nrow(V)) {
    stop("The mean matrix M has different number of parameters than scale
         parameter matrix V: M  -- ", dim(M)[2], "; V -- ", dim(V)[1], ".")
  }
  if (!is.positive.definite(U, tol)) {
    stop("U is not positive definite. Calculation may not be accurate.
         Possibly raise tolerance.")
  }
  if (!is.positive.definite(V, tol)) {
    stop("V is not positive definite. Calculation may not be accurate.
         Possibly raise tolerance.")
  }
  return(invisible())
}

# Unsure if should add
# if (!is.matrix(M)) {  M <- matrix(M) }
