# Convert log abundance covariance matrix to clr covariance matrix
omegaToGamma <- function(omega) {
  p <- ncol(omega)
  G <- diag(p) - matrix(1, p, p)/p
  G %*% omega %*% G
}

# Compute log abundance covariance matrix with given diagonal entries
# which has the given clr covariance matrix
gammaToOmega <- function(gamma, omegaDiag) {
  d <- omegaDiag - diag(gamma)
  R <- matrix(d, length(d), length(d))
  C <- t(R)
  omega <- gamma + 0.5*(R + C)
  if (min(eigen(omega)$values) < 0)
    message("Omega is not a valid covariance matrix")
  omega
}