# Example of alternative covariances and the corresponding
# partial correlations for p = 64

library(corpcor)
source("../cov.R")

Omega <- solve(as.matrix(read.table("../sim/icovs/cluster_64_neg.csv",
                                    sep = ",", header = FALSE)))
p <- ncol(Omega)
Gamma <- omegaToGamma(Omega)
groupLims <- posDefRange(Omega, 1:16)
indivLims <- sapply(1:p, function(i) posDefRange(Omega, i))

par(mfrow = c(2,3))
plotCor(cor2pcor(gammaToOmega(Gamma, diag(Gamma) + 1e-4), "(a)"))
plotCor(cor2pcor(gammaToOmega(Gamma, diag(Gamma) + 1e-3), "(b)"))
plotCor(cor2pcor(gammaToOmega(Gamma, diag(Gamma) + 1e-2), "(c)"))
plotCor(cor2pcor(gammaToOmega(Gamma, c(diag(Omega)[1:16] + 0.99*groupLims[1],
                                       diag(Omega)[17:p])), "(d)"))
plotCor(cor2pcor(gammaToOmega(Gamma, diag(Omega)), "(e)"))
plotCor(cor2pcor(gammaToOmega(Gamma, c(diag(Omega)[1:16] + 0.99*groupLims[2],
                                       diag(Omega)[17:p])), "(f)"))
par(mfrow = c(1,1))
