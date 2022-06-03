library(data.table)
library(scales)

## LOAD AND MERGE PBS DATA

f <- "/home/genis/warthogs/analyses/fst/results3/pbs/desert_southeast_ghana_slidingwindow_100000_stepsize_100000.pbs"

res <- data.table::fread(f, data.table=FALSE, header=TRUE)
#resN <- resN[resN$PBS2>0,]

res$pos <- paste(res$chr, res$midPos, sep="_")

res$start <- res$midPos - 5e4 + 1

res$end <- res$midPos + 5e4

pbs <- res[order(as.integer(gsub("chr","",res$chr)), res$midPos),] # reorder so chrom order is chr1, chr2... rather than chr1, crh10... and then order by coordinate


## LOAD AND MERGE FD DATA
f <- "/home/genis/warthogs/analyses/genome_scan/fd_basedv231052021/abbababa_win100kb_GhanaSEDesert_combined.csv.gz"

fd <- fread(f, sep=",", header=T, data.table=FALSE)

fd <- fd[,-1]
colnames(fd) <- c("chr", colnames(fd)[-1])

res <- merge(x=fd, y=pbs, by=c("chr", "start", "end"), all=F)
res <- res[order(as.integer(gsub("chr","",res$chr)), res$midPos),]
keep <- res$Nsites > quantile(res$Nsites,0.05)
fd <- res

fd_col <- "#F95700FF"
keep2 <- keep & fd$fd > 0


## load list of 19 top hits
tops <- read.table("/home/genis/warthogs/analyses/genome_scan/finalreally_v4_11082021/tophits_fd_100kb_top20RegionsGenis.csv", h=T, sep="\t", stringsAsFactors=F)[1:19,]
ord <- order(as.integer(gsub("chr", "",tops$chr)), tops$midPos)
tops <- tops[ord,]
tops$xInPlot <- apply(cbind(tops$chr, tops$midPos), 1,function(x) which(fd[keep2,]$chr == x[1] & fd[keep2,]$midPos == x[2]))

tops$Label.in.plot <- gsub("MHC complex", "MHC locus", tops$Label.in.plot)


# do not include regions without genes, keep only one MHC complex
k <- !tops$Label.in.plot == "none" & !1:19 %in% c(9, 10)


#outpng <-  "/newHome/genis/mydata/warthogs/analyses/genome_scan/finalreally_v4_11082021/manhattanPlotv2.png"
cex_lab <- 1
#bitmap(outpng, w=9,h=3,res=300)
#par(mar = c(5.1, 4.6, 4.1, 2.1), fig=c(0,1,0,1))
par(mar = c(3.1, 4.6, 2.1, 2.1))
ccol <- c("lightblue","darkblue")

plot(fd$fd[keep2], col= ccol[as.integer(gsub("chr","",fd$chr[keep2]))%%2+1], pch=16, ylab=expression("f"["d"]), xlab="Chromosome", cex.lab=cex_lab, xaxt="n", cex.main=1, bty="n", xpd=NA)

chr <- factor(fd$chr, levels = paste0("chr", 1:18))
tabC <- table(chr[keep2])
mm <- cumsum(tabC)-tabC/2
mtext(1:length(tabC),1,1:length(tabC)%%2,at=unlist(mm), cex=0.6)

abline(h=quantile(fd$fd[keep], 0.999, na.rm=T), lty=2, col=fd_col)

adx <- c(0, -50, 530, 470, 400, 400, 0, 100, 0, 450, 450, 0,0,-100) * 1.3
ady <- c(0, 0, -0.01, -0.04, -0.002, 0, 0, 0, 0, -0.025, -0.01, rep(0, 3)) * 1.1
text(x=tops$xInPlot[k] + adx, cex=1,  y=tops$fd[k]+0.02 + ady, labels=tops$Label.in.plot[k], xpd=NA)
#text(x=tops$xInPlot[k], y=tops$fd[k]+0.02, labels=tops$Label.in.plot[k], xpd=NA)

legend(x=nrow(fd[keep2,])*0.8, y=0.45,bty="n", lty=2, col=fd_col, legend=expression("99.9% quantile f"["d"]), xpd=NA, cex=1)

#par(fig = c(0.84, 1, 0.5, 1), new = T)
hist(fd$fd[keep2], breaks=100, main="", xlab=expression("f"["d"]), ylab="Count", cex.lab=1, xpd=NA)
abline(v=quantile(fd$fd[keep2], 0.999, na.rm=T), lty=2, col=fd_col)
#par(fig = c(0,1,0,1))
#dev.off()
