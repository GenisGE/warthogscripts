source("/home/genis/warthogs/info/loadPopInfo.R")



het_file<- "/home/genis/warthogs/analyses/heterozygosity/results/sfs/collected.txt"

het <- read.table(het_file, h=T)


dat <- merge(df, het, by.x=c("Serial_number"), by.y=c("id"))

dat$Population <- factor(dat$Population, levels=popord)

plotname <- "warthogsHets.png"
bitmap(plotname, width=6, height=6, res=300)
par(oma=c(2,3,0,0))
boxplot(dat$het ~ dat$Population, col=wartCols[popord],outline=F,las=2, xlab="", cex.axis=1.5, cex.lab=1.7)
title(ylab="Heterozygosity", line=5.5, xpd=NA, cex.lab=1.7)
points(dat$het ~ jitter(as.integer(dat$Population),1),pch=21,bg=wartCols[dat$Population])
dev.off()




genoHetsFold <- "/home/genis/warthogs/analyses/psmc/results/hets"
f <- list.files(genoHetsFold, full.names=T)


inds <-  sapply(strsplit(basename(f), "_"), function(x) x[1])
ad <- sapply(strsplit(basename(f), "_"), function(x) gsub(".het","",x[4]))
pops <- c("Zimbabwe", "Namibia", "Zambia", "Ghana", "Tanzania", "Desert", "Ghana", "Ghana")
pops <- do.call("c",lapply(pops,rep,times=3))
pops<- factor(pops, levels=popord)
hets <- sapply(f, scan)

dat2 <- data.frame(id=inds, pop=pops, ad=ad, het=hets, row.names=NULL)
dat2 <- dat2[-grep("group", dat2$id),]

ad_use <- 1
if(FALSE){
plotname <- "warthogsHets2.png"
bitmap(plotname, width=6.5, height=6, res=300)
par(oma=c(0,1,0,0))
boxplot(dat$het ~ dat$Population, col=wartCols[popord],outline=F,las=1, xlab="Population")
title(ylab="Heterozygosity", line=4, xpd=NA)
points(dat$het ~ jitter(as.integer(dat$Population),1),pch=21,bg=wartCols[dat$Population])
points(dat2$het[dat2$ad==ad_use] ~ jitter(as.integer(dat2$pop[dat2$ad==ad_use]),1), pch=4, cex=2, lwd=3, col=wartCols[dat2$pop[dat2$ad==ad_use]])
dev.off()
}




plotname <- "warthogsHetsAllInfo.png"

bitmap(plotname, width=6.5, height=6, res=300)
par(oma=c(2,3,0,0))

boxplot(dat$het ~ dat$Population, col=wartCols[popord],outline=F,las=2, xlab="", cex.axis=1.5, cex.lab=1.7,  ylim=c(0.0015, 0.004))
title(ylab="Heterozygosity", line=5.5, xpd=NA, cex.lab=1.7)

points(dat$het[is.na(dat$High_depth)]~ jitter(as.integer(dat$Population[is.na(dat$High_depth)]),1),pch=21,bg=wartCols[dat$Population[is.na(dat$High_depth)]])
points(dat$het[!is.na(dat$High_depth)]~ jitter(as.integer(dat$Population[!is.na(dat$High_depth)]),1),pch=25,bg=wartCols[dat$Population[!is.na(dat$High_depth)]], cex=2)


points(dat2$het[dat2$ad==1] ~ jitter(as.integer(dat2$pop[dat2$ad==1]),1), pch=3, cex=2, lwd=3, col=wartCols[dat2$pop[dat2$ad==1]])
points(dat2$het[dat2$ad==2] ~ jitter(as.integer(dat2$pop[dat2$ad==2]),1), pch=4, cex=2, lwd=3, col=wartCols[dat2$pop[dat2$ad==2]])
points(dat2$het[dat2$ad==3] ~ jitter(as.integer(dat2$pop[dat2$ad==3]),1), pch=8, cex=2, lwd=3, col=wartCols[dat2$pop[dat2$ad==3]])

legend("topright", legend=c("GL", "CG minad 1", "CG minad 2", "CG minad 3"), pch=c(25, 3,4,8), title="Method for high depth", cex=1.5)

dev.off()

