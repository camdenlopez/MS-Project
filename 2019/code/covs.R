# Generates log(x) covariances from the clr(x)
# covariances of the HMP datasets and computes the
# resulting marginal and partial correlations.

library(plotly)
library(tidyverse)
source("ops.R")
set.seed(547)

load("../data/hmp.RData")

Y.sm <- hmp.sm %>%
  mutate(Count = Count + 1) %>%
  spread(key = OTU, value = Count) %>%
  column_to_rownames(var = "SampleID") %>%
  as.matrix
Y.sm <- Y.sm[, otus$OTU[1:ncol(Y.sm)]]
cov.obs.sm <- cov(log(Y.sm))
cov.clr.sm <- cov_to_clr(cov.obs.sm)
var.log.sm <- generate_covs(cov.clr = cov.clr.sm,
                            var.max = 1.5 * sum(diag(cov.obs.sm)),
                            n.cov = 1e4)
cor.sm <- get_cor(cov.clr.sm, var.log.sm)
pcor.sm <- get_cor(cov.clr.sm, var.log.sm, partial = TRUE)
plot_ly(data.frame(var.log.sm),
        x = ~ var1, y = ~ var2, z = ~ var3) %>%
  add_markers(size = 1)

Y.md <- hmp.md %>%
  mutate(Count = Count + 1) %>%
  spread(key = OTU, value = Count) %>%
  column_to_rownames(var = "SampleID") %>%
  as.matrix
Y.md <- Y.md[, otus$OTU[1:ncol(Y.md)]]
cov.obs.md <- cov(log(Y.md))
cov.clr.md <- cov_to_clr(cov.obs.md)
var.log.md <- generate_covs(cov.clr.md,
                            var.max = 1.5 * sum(diag(cov.obs.md)),
                            n.cov = 1e4)
cor.md <- get_cor(cov.clr.md, var.log.md)
pcor.md <- get_cor(cov.clr.md, var.log.md, partial = TRUE)
plot_ly(data.frame(var.log.md),
        x = ~ var1, y = ~ var2, z = ~ var3) %>%
  add_markers(size = 1)

Y.lg <- hmp.lg %>%
  mutate(Count = Count + 1) %>%
  spread(key = OTU, value = Count) %>%
  column_to_rownames(var = "SampleID") %>%
  as.matrix
Y.lg <- Y.lg[, otus$OTU[1:ncol(Y.lg)]]
cov.obs.lg <- cov(log(Y.lg))
cov.clr.lg <- cov_to_clr(cov.obs.lg)
var.log.lg <- generate_covs(cov.clr.lg,
                            var.max = 1.5 * sum(diag(cov.obs.lg)),
                            n.cov = 1e4)
cor.lg <- get_cor(cov.clr.lg, var.log.lg)
pcor.lg <- get_cor(cov.clr.lg, var.log.lg, partial = TRUE)
plot_ly(data.frame(var.log.lg),
        x = ~ var1, y = ~ var2, z = ~ var3) %>%
  add_markers(size = 1)

save(Y.sm, cov.obs.sm, cov.clr.sm, var.log.sm, cor.sm, pcor.sm,
     Y.md, cov.obs.md, cov.clr.md, var.log.md, cor.md, pcor.md,
     Y.lg, cov.obs.lg, cov.clr.lg, var.log.lg, cor.lg, pcor.lg,
     file = "../data/covs.RData")
