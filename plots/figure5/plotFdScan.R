
outpng <- "fdfigure.png"
png(outpng, h=6, w=8, res=300, units="in")
#bitmap(outpng, h=6, w=8, res=300)
par(oma=c(6,0.5,0,0))
layout(matrix(c(1,1,1,1,2,3,3,5,5,5,4,4,6,6,6), nrow=3, byrow=T), heights=c(0.5,0.2,0.3))
#layout(matrix(c(1,1,1,2,4,4, 3,5,5), nrow=3, byrow=T), heights=c(0.5,0.25,0.25))
source("plotManhattan.R")
source("plotLocusMHCMuc19.R")
dev.off()
