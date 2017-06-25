library(mvtnorm)
library(huge)
library(glasso)
library(igraph)
library(ggplot2)
source("../pr.R")
source("../cov.R")

# Example of potential log abundance covariances

set.seed(656)
W <- as.matrix(read.csv("../../data/otu-table.csv"))
p <- 24
keep <- sample(1:ncol(W), p)
W <- W[,keep]
W <- W[which(rowSums(W) > 0),]
X <- t(apply(W, 1, function(w) w / sum(w)))
Y <- t(apply(X, 1, function(x) x * median(rowSums(W)))) + 1
Z <- t(apply(log(Y), 1, function(x) x - mean(x)))
Gamma <- cov(Z)
Omega <- gammaToOmega(Gamma, diag(Gamma) + 1e-1)
Alt1 <- gammaToOmega(Gamma, diag(Gamma) + 1e-1 + rnorm(p, sd = 0.07))
Alt2 <- gammaToOmega(Gamma, diag(Gamma) + 1e-1 + rnorm(p, sd = 0.07))
Alt3 <- gammaToOmega(Gamma, diag(Gamma) + 1e-1 + rnorm(p, sd = 0.07))

clr.est <- glasso(Gamma, rho = 0.43)
sum(upper.tri(clr.est$wi) & clr.est$wi != 0)
adj <- ifelse(abs(clr.est$wi) > .Machine$double.eps, 1, 0)
g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
lay <- layout_in_circle(g)
E(g)$color <- sapply(1:24, function(i) { 
  e <- get.edgelist(g)[i,]
  ifelse(clr.est$wi[e[1],e[2]] < 0, "blue", "red")
})
pdf("example-clr-cor.pdf", width = 6, height = 6)
par(mar = rep(2, 4))
plotCor(cov2cor(Gamma))
dev.off()
pdf("example-clr-graph.pdf", width = 6, height = 6)
par(mar = rep(2, 4))
plot(g, vertex.size = 8, vertex.color = "black",
     vertex.label = "", edge.width = 3, layout = lay)
dev.off()

alt1.est <- glasso(Alt1, rho = 0.478)
sum(upper.tri(alt1.est$wi) & alt1.est$wi != 0)
adj1 <- ifelse(abs(alt1.est$wi) > .Machine$double.eps, 1, 0)
g1 <- graph_from_adjacency_matrix(adj1, mode = "undirected", diag = FALSE)
E(g1)$color <- sapply(1:24, function(i) { 
  e <- get.edgelist(g)[i,]
  ifelse(alt1.est$wi[e[1],e[2]] < 0, "blue", "red")
})
pdf("example-alt1-cor.pdf", width = 6, height = 6)
par(mar = rep(2, 4))
plotCor(cov2cor(Alt1))
dev.off()
pdf("example-alt1-graph.pdf", width = 6, height = 6)
par(mar = rep(2, 4))
plot(g1, vertex.size = 8, vertex.color = "black",
     vertex.label = "", edge.width = 3, layout = lay)
dev.off()

alt2.est <- glasso(Alt2, rho = 0.48)
sum(upper.tri(alt2.est$wi) & alt2.est$wi != 0)
adj2 <- ifelse(abs(alt2.est$wi) > .Machine$double.eps, 1, 0)
g2 <- graph_from_adjacency_matrix(adj2, mode = "undirected", diag = FALSE)
E(g2)$color <- sapply(1:24, function(i) { 
  e <- get.edgelist(g)[i,]
  ifelse(alt2.est$wi[e[1],e[2]] < 0, "blue", "red")
})
pdf("example-alt2-cor.pdf", width = 6, height = 6)
par(mar = rep(2, 4))
plotCor(cov2cor(Alt2))
dev.off()
pdf("example-alt2-graph.pdf", width = 6, height = 6)
par(mar = rep(2, 4))
plot(g2, vertex.size = 8, vertex.color = "black",
     vertex.label = "", edge.width = 3, layout = lay)
dev.off()

alt3.est <- glasso(Alt3, rho = 0.478)
sum(upper.tri(alt3.est$wi) & alt3.est$wi != 0)
adj3 <- ifelse(abs(alt3.est$wi) > .Machine$double.eps, 1, 0)
g3 <- graph_from_adjacency_matrix(adj3, mode = "undirected", diag = FALSE)
E(g3)$color <- sapply(1:24, function(i) { 
  e <- get.edgelist(g)[i,]
  ifelse(alt3.est$wi[e[1],e[2]] < 0, "blue", "red")
})
pdf("example-alt3-cor.pdf", width = 6, height = 6)
par(mar = rep(2, 4))
plotCor(cov2cor(Alt3))
dev.off()
pdf("example-alt3-graph.pdf", width = 6, height = 6)
par(mar = rep(2, 4))
plot(g3, vertex.size = 8, vertex.color = "black",
     vertex.label = "", edge.width = 3, layout = lay)
dev.off()

# Example of graph selection performance
# with large compositional effect

set.seed(621)
n <- 1024
lam <- seq(from = 0.65, to = 0.001, by = -0.001)

g.cluster <- as.matrix(
  read.table("../sim/graphs/cluster_64.csv",
             sep=",", header=FALSE))
icov.cluster <- as.matrix(
  read.table("../sim/icovs/cluster_64_neg.csv",
             sep=",", header=FALSE))
icov.cluster[1:16,1:16] <- (1/25) * icov.cluster[1:16,1:16]
cov.cluster <- solve(icov.cluster)
X.cluster <- rmvnorm(n, sigma = cov.cluster)
Y.cluster <- t(apply(X.cluster, 1, function(x) x - mean(x)))
X.est.cluster <- huge(X.cluster, method = "glasso", lambda = lam)
Y.est.cluster <- huge(Y.cluster, method = "glasso", lambda = lam)
X.pr.cluster <- precisionRecall(X.est.cluster, g.cluster)
Y.pr.cluster <- precisionRecall(Y.est.cluster, g.cluster)

X.idx <- X.pr.cluster[3,min(which(X.pr.cluster[2,] <= 0.9))]
X.adj <- as.matrix(X.est.cluster$path[[X.idx]])
X.g <- graph_from_adjacency_matrix(X.adj, mode = "undirected", diag = FALSE)
lay <- layout_in_circle(X.g)
E(X.g)$color <- sapply(1:sum(X.adj[upper.tri(X.adj)]), function(i) { 
  e <- get.edgelist(X.g)[i,]
  ifelse(X.est.cluster$icov[[X.idx]][e[1],e[2]] < 0, "blue", "red")
})

e <- sum(X.adj[upper.tri(X.adj)])
Y.idx <- min(which(sapply(Y.est.cluster$path, function(g) sum(g) / 2) >= e))

Y.adj <- as.matrix(Y.est.cluster$path[[Y.idx]])
Y.g <- graph_from_adjacency_matrix(Y.adj, mode = "undirected", diag = FALSE)
E(Y.g)$color <- sapply(1:sum(Y.adj[upper.tri(Y.adj)]), function(i) { 
  e <- get.edgelist(Y.g)[i,]
  ifelse(Y.est.cluster$icov[[Y.idx]][e[1],e[2]] < 0, "blue", "red")
})

true.g <- graph_from_adjacency_matrix(g.cluster, mode = "undirected", diag = FALSE)
E(true.g)$color <- "blue"
pdf("cluster-graph.pdf", width = 6, height = 6)
par(mar = rep(2,4))
plot(true.g, vertex.size = 6, vertex.color = c(rep("yellow", 16), rep("black", 48)),
     vertex.label = "", edge.width = 3, layout = lay)
dev.off()
pdf("cluster-log.pdf", width = 6, height = 6)
par(mar = rep(2,4))
plot(X.g, vertex.size = 6, vertex.color = c(rep("yellow", 16), rep("black", 48)),
     vertex.label = "", edge.width = 3, layout = lay)
dev.off()
pdf("cluster-clr.pdf", width = 6, height = 6)
par(mar = rep(2,4))
plot(Y.g, vertex.size = 6, vertex.color = c(rep("yellow", 16), rep("black", 48)),
     vertex.label = "", edge.width = 3, layout = lay)
dev.off()

pr <- data.frame(Recall = c(X.pr.cluster[1,], Y.pr.cluster[1,]),
                 Precision = c(X.pr.cluster[2,], Y.pr.cluster[2,]),
                 Data = factor(c(rep("log", ncol(X.pr.cluster)),
                                 rep("clr", ncol(Y.pr.cluster)))))
pr <- pr[pr$Precision >= 0.04 | pr$Recall == 1,]
pdf("cluster-pr.pdf", width = 6, height = 4)
par(mar = rep(2,4))
ggplot(pr, aes(x=Recall, y=Precision, linetype=Data)) +
  geom_line() + theme_minimal(base_size = 16)
dev.off()