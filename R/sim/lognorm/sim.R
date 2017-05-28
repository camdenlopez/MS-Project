library(mvtnorm)
library(huge)
library(foreach)
library(doMC)

source("../../pr.R")
source("../../cov.R")

registerDoMC(cores = 4)
set.seed(940)

reps <- 20
method.levels <- c("glasso", "mb")
p.levels <- c(64, 256)
graph.levels <- c("band", "cluster", "scalefree")
icov.levels <- c("neg", "mix", "pos")
n.levels <- c(32, 64, 128, 256, 512)
n.sim <- length(method.levels) * length(p.levels) *
  length(graph.levels) * length(icov.levels) *
  length(n.levels) * reps
design <- data.frame(
  method = gl(length(method.levels),
              length(p.levels) * length(graph.levels) *
                length(icov.levels) * length(n.levels) * reps,
              n.sim),
  p = gl(length(p.levels),
         length(graph.levels) * length(icov.levels) *
           length(n.levels) * reps,
         n.sim),
  graph = gl(length(graph.levels),
             length(icov.levels) * length(n.levels) * reps,
             n.sim),
  icov = gl(length(icov.levels),
            length(n.levels) * reps,
            n.sim),
  n = gl(length(n.levels),
         reps,
         n.sim))
levels(design$method) <- method.levels
levels(design$p) <- p.levels
levels(design$graph) <- graph.levels
levels(design$icov) <- icov.levels
levels(design$n) <- n.levels

start <- 1
end <- n.sim
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  start <- as.integer(args[1])
  end <- as.integer(args[2])
}

metrics <- c("log.recall.min", "log.recall.max", "log.aupr",
             "clr.recall.min", "clr.recall.max", "clr.aupr")

results <- foreach(i = start:end, .combine = rbind) %dopar% {
  method <- method.levels[design$method[i]]
  p <- p.levels[design$p[i]]
  graph <- graph.levels[design$graph[i]]
  icov <- icov.levels[design$icov[i]]
  n <- n.levels[design$n[i]]
  
  message(sprintf("Trial %d: method=%s p=%d graph=%s icov=%s n=%d",
                  i, method, p, graph, icov, n))
  
  trial <- numeric(length(metrics))
  names(trial) <- metrics
  
  g <- as.matrix(read.table(sprintf("../graphs/%s_%d.csv",
                                    graph, p),
                            sep=",", header=FALSE))
  prec <- as.matrix(read.table(sprintf("../icovs/%s_%d_%s.csv",
                                       graph, p, icov),
                               sep=",", header=FALSE))
  cov <- solve(prec)
  
  X <- mvtnorm::rmvnorm(n, sigma = cov)
  log.est <- huge(X, method = method, nlambda = 30,
                  lambda.min.ratio = 1/100, verbose = FALSE)
  
  log.pr <- precisionRecall(log.est, g)
  trial["log.recall.min"] <- min(log.pr[1,])
  trial["log.recall.max"] <- max(log.pr[1,])
  trial["log.aupr"] <- areaUnderPR(log.pr)
  
  Y <- t(apply(X, 1, function(x) x - mean(x)))
  clr.est <- huge(Y, method = "glasso", nlambda = 30,
                  lambda.min.ratio = 1/100, verbose = FALSE)
  
  clr.pr <- precisionRecall(clr.est, g)
  trial["clr.recall.min"] <- min(clr.pr[1,])
  trial["clr.recall.max"] <- max(clr.pr[1,])
  trial["clr.aupr"] <- areaUnderPR(clr.pr)
  
  trial
}

results <- cbind(design[start:end,], results)
write.csv(results, "results.csv", row.names = FALSE)