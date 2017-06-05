# Convert log abundance covariance matrix to clr covariance matrix
omegaToGamma <- function(omega) {
  p <- ncol(omega)
  G <- diag(p) - matrix(1, p, p)/p
  G %*% omega %*% G
}

# Compute log abundance covariance matrix with given diagonal entries
# and given clr covariance matrix
gammaToOmega <- function(gamma, omegaDiag) {
  d <- omegaDiag - diag(gamma)
  R <- matrix(d, length(d), length(d))
  C <- t(R)
  gamma + 0.5*(R + C)
}

# Find range of diagonal entry values for which the potential
# log abundance covariance matrix is positive definite, varying
# only those entries specified by idx while holding others fixed.
posDefRange <- function(omega, idx, step = 1e-1, tol = 1e-4) {
  Gamma <- omegaToGamma(omega)
  # LL = lower boundary of lower limit for valid entries
  # LU = upper boundary of lower limit for valid entries
  LL <- LU <- 0
  # First find lower boundary using exponentially
  # growing step size
  while (TRUE) {
    D <- diag(omega)
    D[idx] <- D[idx] + LL
    Omega <- gammaToOmega(Gamma, D)
    if (min(eigen(Omega)$values) < 0)
      break
    else {
      LU <- LL
      LL <- LL - step
      step <- 2 * step
    }
  }
  # Do binary search to find lower limit within the
  # given tolerance
  while (abs(LU - LL) > tol) {
    D <- diag(omega)
    D[idx] <- D[idx] + (LL + LU) / 2
    Omega <- gammaToOmega(Gamma, D)
    if (min(eigen(Omega)$values) < 0)
      LL <- (LL + LU) / 2
    else
      LU <- (LL + LU) / 2
  }
  # Do same procedure as above, except for upper limit
  UL <- UU <- 0
  while (TRUE) {
    D <- diag(omega)
    D[idx] <- D[idx] + UU
    Omega <- gammaToOmega(Gamma, D)
    if (min(eigen(Omega)$values) < 0)
      break
    else {
      UL <- UU
      UU <- UU + step
      step <- 2 * step
    }
  }
  while (abs(UU - UL) > tol) {
    D <- diag(omega)
    D[idx] <- D[idx] + (UL + UU) / 2
    Omega <- gammaToOmega(Gamma, D)
    if (min(eigen(Omega)$values) < 0)
      UU <- (UL + UU) / 2
    else
      UL <- (UL + UU) / 2
  }
  c(LU, UL)
}