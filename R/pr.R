# Compute precision-recall curve from solution path returned by "huge"
precisionRecall <- function(est, graph) {
  pr <- sapply(1:length(est$lambda), function(i) {
    E <- as.matrix(est$path[[i]])[upper.tri(as.matrix(est$path[[i]]))]
    G <- as.matrix(graph)[upper.tri(as.matrix(graph))]
    sum(E * G) / c(sum(G), sum(E))
  })
  pr <- pr[,!is.na(pr[1,]) & !is.na(pr[2,])]
}

# Compute area under precision-recall curve
areaUnderPR <- function(pr) {
  if (pr[1,1] > 0)
    pr <- cbind(c(0, pr[2,1]), pr)
  if (pr[1,ncol(pr)] < 1)
    pr <- cbind(pr, c(1, pr[2,ncol(pr)]))
  sum(diff(pr[1,]) * (pr[2,-ncol(pr)] + pr[2,-1]) / 2)
}