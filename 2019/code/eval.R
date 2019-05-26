library(tidyverse)
theme_set(theme_light())

load("../data/covs.RData")

cor.min.sm <- apply(cor.sm, 2, min)
cor.max.sm <- apply(cor.sm, 2, max)
table(cor.min.sm < 0 & cor.max.sm > 0)

pcor.min.sm <- apply(pcor.sm, 2, min)
pcor.max.sm <- apply(pcor.sm, 2, max)
table(pcor.min.sm < 0 & pcor.max.sm > 0)

cor.min.md <- apply(cor.md, 2, min)
cor.max.md <- apply(cor.md, 2, max)
table(cor.min.md < 0 & cor.max.md > 0)

pcor.min.md <- apply(pcor.md, 2, min)
pcor.max.md <- apply(pcor.md, 2, max)
table(pcor.min.md < 0 & pcor.max.md > 0)

cor.min.lg <- apply(cor.lg, 2, min)
cor.max.lg <- apply(cor.lg, 2, max)
table(cor.min.lg < 0 & cor.max.lg > 0)

pcor.min.lg <- apply(pcor.lg, 2, min)
pcor.max.lg <- apply(pcor.lg, 2, max)
table(pcor.min.lg < 0 & pcor.max.lg > 0)
