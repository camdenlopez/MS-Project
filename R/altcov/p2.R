source("../cov.R")

omega22lims <- function(omega11, gamma11) {
  omega11 + 4*gamma11 + c(-1,1)*4*sqrt(gamma11*omega11)
}

plotOmegaRegion <- function(gamma11, lims, shadeAlpha = 1) {
  omega11 <- seq(lims[1], lims[2], length.out = 1000)
  omega22 <- sapply(omega11, function(x) omega22lims(x, gamma11))
  omega22 <- apply(omega22, 2, function(x) c(x[1], min(c(x[2], lims[2]))))
  plot(omega11, omega22[1,], type = "l",
       xlim = lims, ylim = lims,
       xlab = expression(omega[11]), ylab = expression(omega[22]),
       main = substitute(paste(gamma[11], " = ", g11, sep = ""),
                         list(g11 = gamma11)))
  lines(omega11[omega22[2,] < lims[2]], omega22[2,omega22[2,] < lims[2]])
  polygon(c(rev(omega11), omega11), c(rev(omega22[1,]), omega22[2,]),
          density = 15, angle = -15, lty = 1, border = NA,
          col = rgb(0, 0, 0, shadeAlpha))
}

pdf("p2-regions.pdf", width=16, height=6)
par(mfrow=c(1,3), cex=1.6)
plotOmegaRegion(0.1, c(0, 30))
plotOmegaRegion(  1, c(0, 30))
plotOmegaRegion(  1, c(0,  5), shadeAlpha=0.2)
lines(c(0,4), c(4,0), lty = 2)
text(x = 1.5*gamma11, y = 1.5*gamma11, labels = expression(omega[12]<0))
text(x = 2.75*gamma11, y = 2.75*gamma11, labels = expression(omega[12]>0))
par(mfrow = c(1,1), cex=1)
dev.off()