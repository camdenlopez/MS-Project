library(SpiecEasi)
library(mvtnorm)
library(huge)

source("../../pr.R")

set.seed(356)

g.band <- as.matrix(
  read.table("../graphs/band_64.csv",
             sep=",", header=FALSE))
g.cluster <- as.matrix(
  read.table("../graphs/cluster_64.csv",
             sep=",", header=FALSE))
g.scalefree <- as.matrix(
  read.table("../graphs/scalefree_64.csv",
             sep=",", header=FALSE))
cov.band <- solve(
  read.table("../icovs/band_64_mix.csv",
             sep=",", header=FALSE))
cov.cluster <- solve(
  read.table("../icovs/cluster_64_mix.csv",
             sep=",", header=FALSE))
cov.scalefree <- solve(
  read.table("../icovs/scalefree_64_mix.csv",
             sep=",", header=FALSE))

M <- 400
n <- 256

cov.band[1,1] <- M * cov.band[1,1]
min(eigen(cov.band)$values)
X.band <- rmvnorm(n, sigma = cov.band)
Y.band <- t(apply(X.band, 1, function(x) x - mean(x)))
X.est.band <- huge(X.band, method = "glasso", nlambda = 30, lambda.min.ratio = 1/100)
Y.est.band <- huge(Y.band, method = "glasso", nlambda = 30, lambda.min.ratio = 1/100)
X.pr.band <- precisionRecall(X.est.band, g.band)
Y.pr.band <- precisionRecall(Y.est.band, g.band)

cov.cluster[1,1] <- M * cov.cluster[1,1]
min(eigen(cov.cluster)$values)
X.cluster <- rmvnorm(n, sigma = cov.cluster)
Y.cluster <- t(apply(X.cluster, 1, function(x) x - mean(x)))
X.est.cluster <- huge(X.cluster, method = "glasso", nlambda = 30, lambda.min.ratio = 1/100)
Y.est.cluster <- huge(Y.cluster, method = "glasso", nlambda = 30, lambda.min.ratio = 1/100)
X.pr.cluster <- precisionRecall(X.est.cluster, g.cluster)
Y.pr.cluster <- precisionRecall(Y.est.cluster, g.cluster)

cov.scalefree[1,1] <- M * cov.scalefree[1,1]
min(eigen(cov.scalefree)$values)
X.scalefree <- rmvnorm(n, sigma = cov.scalefree)
Y.scalefree <- t(apply(X.scalefree, 1, function(x) x - mean(x)))
X.est.scalefree <- huge(X.scalefree, method = "glasso", nlambda = 30, lambda.min.ratio = 1/100)
Y.est.scalefree <- huge(Y.scalefree, method = "glasso", nlambda = 30, lambda.min.ratio = 1/100)
X.pr.scalefree <- precisionRecall(X.est.scalefree, g.scalefree)
Y.pr.scalefree <- precisionRecall(Y.est.scalefree, g.scalefree)

pdf("plots/var-dom-band.pdf", width = 16, height = 6)
par(mfrow = c(1,3), cex = 1.6)
plot(adj2igraph(g.band),
     vertex.size=3, vertex.color="black", vertex.label="",
     edge.color="black", layout=layout_in_circle)
plot(adj2igraph(X.est.band$path[[5]]),
     vertex.size=3, vertex.color="black", vertex.label="",
     edge.color="black", layout=layout_in_circle)
plot(adj2igraph(Y.est.band$path[[5]]),
     vertex.size=3, vertex.color="black", vertex.label="",
     edge.color="black", layout=layout_in_circle)
dev.off()

pdf("plots/var-dom-cluster.pdf", width = 16, height = 6)
par(mfrow = c(1,3), cex = 1.6)
plot(adj2igraph(g.cluster),
     vertex.size=3, vertex.color="black", vertex.label="",
     edge.color="black", layout=layout_in_circle)
plot(adj2igraph(X.est.cluster$path[[5]]),
     vertex.size=3, vertex.color="black", vertex.label="",
     edge.color="black", layout=layout_in_circle)
plot(adj2igraph(Y.est.cluster$path[[5]]),
     vertex.size=3, vertex.color="black", vertex.label="",
     edge.color="black", layout=layout_in_circle)
dev.off()

pdf("plots/var-dom-scalefree.pdf", width = 16, height = 6)
par(mfrow = c(1,3), cex = 1.6)
coords <- layout_(adj2igraph(g.scalefree), nicely())
plot(adj2igraph(g.scalefree),
     vertex.size=3, vertex.color="black", vertex.label="",
     edge.color="black", layout=coords)
plot(adj2igraph(X.est.scalefree$path[[5]]),
     vertex.size=3, vertex.color="black", vertex.label="",
     edge.color="black", layout=coords)
plot(adj2igraph(Y.est.scalefree$path[[5]]),
     vertex.size=3, vertex.color="black", vertex.label="",
     edge.color="black", layout=coords)
dev.off()

pdf("plots/var-dom-pr.pdf", width = 16, height = 6)
par(mfrow = c(1,3), cex = 1.6, lwd = 3)

plot(X.pr.band[1,], X.pr.band[2,], xlim=0:1, ylim=0:1, type="l", lty = 2,
     xlab="Recall", ylab="Precision", axes = FALSE, main = "Band")
axis(1, at=c(0, 0.5, 1), labels=c(0,0.5,1))
axis(2, at=c(0, 0.5, 1), labels=c(0,0.5,1))
points(Y.pr.band[1,], Y.pr.band[2,], type="l", pch=2, lty=1)

plot(X.pr.cluster[1,], X.pr.cluster[2,], xlim=0:1, ylim=0:1, type="l", lty = 2,
     xlab="Recall", ylab="Precision", axes = FALSE, main = "Cluster")
axis(1, at=c(0, 0.5, 1), labels=c(0,0.5,1))
axis(2, at=c(0, 0.5, 1), labels=c(0,0.5,1))
points(Y.pr.cluster[1,], Y.pr.cluster[2,], type="l", pch=2, lty=1)

plot(X.pr.scalefree[1,], X.pr.scalefree[2,], xlim=0:1, ylim=0:1, type="l", lty = 2,
     xlab="Recall", ylab="Precision", axes = FALSE, main = "Scale-Free")
axis(1, at=c(0, 0.5, 1), labels=c(0,0.5,1))
axis(2, at=c(0, 0.5, 1), labels=c(0,0.5,1))
points(Y.pr.scalefree[1,], Y.pr.scalefree[2,], type="l", pch=2, lty=1)

dev.off()