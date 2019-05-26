# Functions implementing various operations
# and algorithms for covariance matrices.

library(corpcor)
library(corrplot)

# Convert log(x) covariance matrix cov.log to the
# corresponding clr(x) covariance matrix.
cov_to_clr <- function (cov.log) {
  D <- nrow(cov.log)
  G <- diag(D) - (1 / D)
  G %*% cov.log %*% G
}

# Convert clr(x) covariance matrix cov.clr to the
# corresponding log(x) covariance matrix that has
# variances var.log.
cov_to_log <- function (cov.clr, var.log) {
  D <- nrow(cov.clr)
  var.diff <- matrix(rep(var.log - diag(cov.clr), D),
                     nrow = D, byrow = TRUE)
  cov.log <- cov.clr + 0.5 * (var.diff + t(var.diff))
  if (is.positive.definite(cov.log)) {
    cov.log
  } else {
    NA
  }
}

# Generate random D-dimensional unit-length vector.
rand_direction <- function (D) {
  dir <- rnorm(D)
  dir / sqrt(sum(dir ^ 2))
}

# Generate n.cov random log(x) covariances which would
# have clr(x) covariances cov.clr, limited to having
# a sum of log(x) variances no greater than var.max.
generate_covs <- function (cov.clr, var.max, n.cov) {
  if (var.max < sum(diag(cov.clr)))
    stop("Cannot generate covs: var.max < tr(cov.clr)")
  D <- nrow(cov.clr)
  var.log <- replicate(n.cov, {
    a <- rand_direction(D)
    A <- matrix(a, nrow = D, ncol = D, byrow = FALSE)
    J <- matrix(1, nrow = D, ncol = D)
    j <- rep(1, D)
    quadform.a <- 0.5 * D * (a %*% cov.clr %*% a)[1, 1]
    quadform.b <- (a %*% cov.clr %*% j)[1, 1]
    quadform.c <- -0.5 * (var.max - sum(diag(cov.clr)))
    r.max <-
      (-quadform.b +
         sqrt(quadform.b ^ 2 - 4 * quadform.a * quadform.c)) /
      (2 * quadform.a)
    r <- runif(1, min = 0, max = r.max)
    cov.log <- cov.clr +
      r * (cov.clr %*% A) +
      r * (t(A) %*% cov.clr) +
      (r ^ 2) * (a %*% cov.clr %*% a)[1, 1] * J
    var.z <- runif(1, min = 0, max = (var.max - sum(diag(cov.log))) / D)
    cov.log <- cov.log + var.z * J
    diag(cov.log)
  }) %>% t
  colnames(var.log) <- paste0("var", 1:D)
  var.log
}

# Plot marginal or partial correlations of
# the log(x) covariance matrix cov.log.
plot_cor <- function (cov.log, partial = FALSE) {
  colnames(cov.log) <- rownames(cov.log) <- NULL
  if (partial) {
    corrplot(cor2pcor(cov.log), method = "color")
  } else {
    corrplot(cov2cor(cov.log), method = "color")
  }
}

# Compute marginal or partial correlations for
# each log(x) covariance matrix corresponding
# to clr(x) covariances cov.clr and log(x)
# variances specified by a row in var.log.
get_cor <- function (cov.clr, var.log, partial = FALSE) {
  ij <- t(combn(1:ncol(var.log), 2))
  cors <- apply(var.log, 1, function (x) {
    cov.log <- cov_to_log(cov.clr, x)
    if (partial) {
      cor2pcor(cov.log)[ij]
    } else {
      cov2cor(cov.log)[ij]
    }
  })
  cors <- data.frame(t(cors))
  names(cors) <- paste(ij[, 1], ij[, 2], sep = ".")
  cors
}
