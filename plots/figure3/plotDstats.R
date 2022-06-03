
# load format d stats based on high depth samples
f <- "/home/genis/warthogs/analyses/dstats_pig/dstat_qpdstat2/dstats.out"

dstat1 <- read.table(f, skip=23, stringsAsFactors=F)[,-1]
colnames(dstat1) <- c("P1", "P2", "P3", "Outgroup", "Dstat", "StdErr", "Z", "nBABA", "nABBA", "nSNPs")
dstat1[, c("Dstat", "Z")] <- dstat1[, c("Dstat", "Z")] * - 1


# load and format dstats based on low depth samples
f <- "/home/genis/warthogs/analyses/dstats_pig/results/dstats/out_dstats.txt"
source("/home/genis/warthogs/info/loadPopInfo.R")
source("/home/genis/warthogs/analyses/dstats_pig/getDstats.R")
library(data.table)

dat <- fread(f, h=T, data.table=F, sep="\t")
dat <- dat[,-ncol(dat)]

pop <- df$Population
names(pop) <- df$Serial_number

popord <- c("Ghana", "Kenya", "Tanzania", "Zimbabwe", "Namibia")
pairs <- t(combn(popord,2))
poptrees <- cbind(pairs, rep("Desert", nrow(pairs)))

hidepth <- df$Serial_number[df$High_depth == "yes" & !is.na(df$High_depth)]


dstats2 <- do.call(rbind,apply(poptrees, 1, getDstat3pop, dstatdf=dat, method=1, locIdx=pop, exclude= hidepth))





outpng <- "/home/genis/warthogs/paper_plots/dstats/dstatsAll.png"
bitmap(outpng, w=8, h=8, res=300)
# now plot both
pops <- c("Ghana", "Kenya", "Tanzania", "Zambia", "Zimbabwe", "Namibia")
pairs <- t(combn(pops,2))
# remove kenya zambia comparison since there is no overlapping data
pairs <- pairs[!(pairs[,1] == "Kenya" & pairs[,2] == "Zambia"),]
plabs <- apply(X = cbind(paste("P1 =", pairs[,1]), paste("P2 =", pairs[,2])), MARGIN=1, FUN=paste, collapse="\n")

# VERTICAL FORMAT X AXIS D STAT Y AXIS POP
par(oma=c(0,3,8,0), xpd=F)
plot(1, type="n", ylab="", xlab="D statistic", xlim=c(-0.1, 0.2), ylim=c(1, nrow(pairs)), yaxt="n", cex.lab=1)
axis(side=2, at=1:nrow(pairs), labels=plabs, las=2, xpd=NA, cex.axis=0.8)
abline(v=0,lty=2)


# get coords to plot low depth, first for boxplot, then for points
at2 <- sapply(gsub("[.]", " ", gsub(".Desert","",unique(dstats2$poptree))), function(x) which(x==apply(pairs,1,paste, collapse=" ")))
at3 <- sapply(gsub("[.]", " ", gsub(".Desert","",dstats2$poptree)), function(x) which(x==apply(pairs,1,paste, collapse=" ")))
dstats2$pair <- factor(gsub("[.]", " ", gsub(".Desert","",dstats2$poptree)), levels=unique(gsub("[.]", " ", gsub(".Desert","",dstats2$poptree))))


boxplot(dstats2$Dstat ~ dstats2$pair, add=T, at=at2, horizontal=T, ylab="", yaxt="n", border="lightgrey", boxwex=0.5)
points(x=dstats2$Dstat, y=jitter(at3, 0.5), col=ifelse(abs(dstats2$Z)>3, "goldenrod", "lightgrey"),
       pch=ifelse(dstats2$H1 == "8931" | dstats2$H2 == "8931", 4, 20), cex=0.75)


# get coords to plot high depth values
at1 <- sapply(paste(dstat1$P1, dstat1$P2), function(x) which(x==apply(pairs,1,paste, collapse=" ")))

#plot high depth d stat vlaues with standard errors
points(y=at1, x=dstat1$Dstat, pch=19, cex=2)
#segments(y0=at1, y1=at1, x0=dstat1$Dstat-dstat1$StdErr, x1=dstat1$Dstat + dstat1$StdErr, lwd=3)
segments(y0=at1, y1=at1, x0=dstat1$Dstat-dstat1$StdErr * 3 , x1=dstat1$Dstat + dstat1$StdErr * 3, lwd=2)

# add legend
leg <- expression(paste("Dstat", phantom(.)%+-%phantom(.), "3 SE"))

legend(x=0.05, y=nrow(pairs), bty="n", pch=19,cex=1,lty=1,lwd=2,legend= leg, title="Genotye calls based")

legend(x=0.05, y=nrow(pairs) - 2.5, legend=c("Dstat (|Z score| < 3)", "Dstat (|Z score| > 3)"), col=c("lightgrey", "goldenrod"), pch=20, xpd=NA,bty='n', cex=1, pt.cex=1.5, title="Single read sampling based")


legend(x=0.05, y=nrow(pairs) - 5, legend = "Kenya sample from range\noverlap with desert warthog", pch=4, xpd=NA,bty='n', cex=1, pt.cex=1.5)

#legend("topright", bty="n", pch=19,col=c("black","goldenrod", "lightgrey"), cex=1,pt.cex=c(2,1,1), lty=c(1, NA,NA),lwd=c(2,NA,NA),legend=c("Dstat +- 3 SD", "Dstat |Z| > 3", "Dstat |Z| < 3")

par(xpd=NA)
lwd=2

# plot right tree
y0 <- nrow(pairs) + 1
y1 <- nrow(pairs) + 4

midx <- 0.1 # middle x postion where root is drawn
w <- 0.05 # distance of mid where tree expands right and left
 
lines(c(midx,midx-w), c(y0,y1), lwd=lwd)
lines(c(midx,midx+w), c(y0,y1), lwd=lwd)

lines(c(midx-w*0.5,midx),c(y0+1.5,y1),lwd=lwd)
lines(c(midx-w * 0.75, midx-w*0.5),c(y0+2.25,y1),lwd=lwd)
text(labels=c("P1","P2","Desert\nwarthog","Domestic\npig"),y=y1+0.5,x=c(midx-w, midx-w*0.5, midx,midx+w),xpd=NA, cex=1)

# add arrow
arrows(x1= midx-w*0.5, y0=y1-0.25, x0= midx - w *0.15, y1=y1-0.25, col='darkred', length=0.1,lwd=2)


# plot tree in left
y0 <- nrow(pairs) + 1
y1 <- nrow(pairs) + 4

midx <- - 0.05 # middle x postion where root is drawn
w <- 0.05 # distance of mid where tree expands right and left
 
lines(c(midx,midx-w), c(y0,y1), lwd=lwd)
lines(c(midx,midx+w), c(y0,y1), lwd=lwd)

lines(c(midx-w*0.5,midx),c(y0+1.5,y1),lwd=lwd)
lines(c(midx-w * 0.75, midx-w*0.5),c(y0+2.25,y1),lwd=lwd)
text(labels=c("P1","P2","Desert\nwarthog","Domestic\npig"),y=y1+0.5,x=c(midx-w, midx-w*0.5, midx,midx+w),xpd=NA, cex=1)

# add arrow
arrows(x1= midx-w*0.85, y0=y1-0.25, x0= midx - w *0.15, y1=y1-0.25, col='darkred', length=0.1,lwd=2)

dev.off()





if(FALSE){
    
bitmap(outpng, w=8, h=8, res=300)
# now plot both
pops <- c("Ghana", "Kenya", "Tanzania", "Zambia", "Zimbabwe", "Namibia")
pairs <- t(combn(pops,2))
# remove kenya zambia comparison since there is no overlapping data
pairs <- pairs[!(pairs[,1] == "Kenya" & pairs[,2] == "Zambia"),]
plabs <- apply(X = cbind(paste("P1 =", pairs[,1]), paste("P2 =", pairs[,2])), MARGIN=1, FUN=paste, collapse="\n")

# VERTICAL FORMAT X AXIS D STAT Y AXIS POP
par(oma=c(0,3,8,0), xpd=F)
plot(1, type="n", ylab="", xlab="D statistic", xlim=c(-0.2, 0.2), ylim=c(1, nrow(pairs)), yaxt="n", cex.lab=1.5)
axis(side=2, at=1:nrow(pairs), labels=plabs, las=2, xpd=NA)
abline(v=0,lty=2)


# get coords to plot low depth, first for boxplot, then for points
at2 <- sapply(gsub("[.]", " ", gsub(".Desert","",unique(dstats2$poptree))), function(x) which(x==apply(pairs,1,paste, collapse=" ")))
at3 <- sapply(gsub("[.]", " ", gsub(".Desert","",dstats2$poptree)), function(x) which(x==apply(pairs,1,paste, collapse=" ")))
dstats2$pair <- factor(gsub("[.]", " ", gsub(".Desert","",dstats2$poptree)), levels=unique(gsub("[.]", " ", gsub(".Desert","",dstats2$poptree))))


boxplot(dstats2$Dstat ~ dstats2$pair, add=T, at=at2, horizontal=T, ylab="", yaxt="n", border="lightgrey", boxwex=0.5)
points(x=dstats2$Dstat, y=jitter(at3, 0.5), col=ifelse(abs(dstats2$Z)>3, "goldenrod", "lightgrey"), pch=20, cex=0.75)


# get coords to plot high depth values
at1 <- sapply(paste(dstat1$P1, dstat1$P2), function(x) which(x==apply(pairs,1,paste, collapse=" ")))

#plot high depth d stat vlaues with standard errors
points(y=at1, x=dstat1$Dstat, pch=19, cex=2)
#segments(y0=at1, y1=at1, x0=dstat1$Dstat-dstat1$StdErr, x1=dstat1$Dstat + dstat1$StdErr, lwd=3)
segments(y0=at1, y1=at1, x0=dstat1$Dstat-dstat1$StdErr * 3 , x1=dstat1$Dstat + dstat1$StdErr * 3, lwd=2)

# add legend
leg <- expression(paste("Dstat", phantom(.)%+-%phantom(.), "3 SE"))

legend(x=0.1, y=nrow(pairs), bty="n", pch=19,cex=1,lty=1,lwd=2,legend= leg, title="Genotye calls based")

legend(x=0.1, y=nrow(pairs)-1.5, legend=c("Dstat (|Z score| < 3)", "Dstat (|Z score| > 3)"), col=c("lightgrey", "goldenrod"), pch=20, xpd=NA,bty='n', cex=1, pt.cex=1, title="Single read sampling based")

#legend("topright", bty="n", pch=19,col=c("black","goldenrod", "lightgrey"), cex=1,pt.cex=c(2,1,1), lty=c(1, NA,NA),lwd=c(2,NA,NA),legend=c("Dstat +- 3 SD", "Dstat |Z| > 3", "Dstat |Z| < 3")


par(xpd=NA)
lwd=2

# plot right tree
y0 <- nrow(pairs) + 1
y1 <- nrow(pairs) + 4

midx <- 0.1 # middle x postion where root is drawn
w <- 0.05 # distance of mid where tree expands right and left
 
lines(c(midx,midx-w), c(y0,y1), lwd=lwd)
lines(c(midx,midx+w), c(y0,y1), lwd=lwd)

lines(c(midx-w*0.5,midx),c(y0+1.5,y1),lwd=lwd)
lines(c(midx-w * 0.75, midx-w*0.5),c(y0+2.25,y1),lwd=lwd)
text(labels=c("P1","P2","Desert\nwarthog","Domestic\npig"),y=y1+0.5,x=c(midx-w, midx-w*0.5, midx,midx+w),xpd=NA)

# add arrow
arrows(x1= midx-w*0.5, y0=y1-0.25, x0= midx - w *0.15, y1=y1-0.25, col='darkred', length=0.1,lwd=2)


# plot tree in left
y0 <- nrow(pairs) + 1
y1 <- nrow(pairs) + 4

midx <- - 0.1 # middle x postion where root is drawn
w <- 0.05 # distance of mid where tree expands right and left
 
lines(c(midx,midx-w), c(y0,y1), lwd=lwd)
lines(c(midx,midx+w), c(y0,y1), lwd=lwd)

lines(c(midx-w*0.5,midx),c(y0+1.5,y1),lwd=lwd)
lines(c(midx-w * 0.75, midx-w*0.5),c(y0+2.25,y1),lwd=lwd)
text(labels=c("P1","P2","Desert\nwarthog","Domestic\npig"),y=y1+0.5,x=c(midx-w, midx-w*0.5, midx,midx+w),xpd=NA)

# add arrow
arrows(x1= midx-w*0.85, y0=y1-0.25, x0= midx - w *0.15, y1=y1-0.25, col='darkred', length=0.1,lwd=2)

dev.off()


}

