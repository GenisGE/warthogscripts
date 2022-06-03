library(data.table)
library(scales)

library(RColorBrewer)

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

fd <- fread(f, sep=",", header=T)

fd <- fd[,-1]
colnames(fd) <- c("chr", colnames(fd)[-1])

res <- merge(x=fd, y=pbs, by=c("chr", "start", "end"), all=F)

res <- res[res$Nsites > quantile(res$Nsites,0.05),]
#res <- res[res$sitesUsed > 100, ]



## to add annotation
dat<- fread("/pontusData/genis/reference/pig_annotation/GCF_000003025.6_Sscrofa11.1_genomic.geneRef",head=F, data.table=F)
colnames(dat) <- c("name", "chr", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
dat$name2 <- gsub("gene-", "", dat$name2)


# load mask bed
f <- "/home/genis/warthogs/genome_masks/all/beds/mask.bed"
mask <- fread(f, h=F, col.names=c("chr", "start", "end"))


q_fst <- quantile(x=res$Fst01, 0.005, na.rm=T)
q_fd <- quantile(x=res$fd, 0.999, na.rm=T)


plotSegment <- function(chr, midpos, winplot=2e6){

    c <- chr
    poss <- c(midpos-winplot/2, midpos+winplot/2)

    par(mar=c(0, 4.2, 2.1, 2.1))    

    xforrs = 0.03
    regsz = 1.2
    width=22 
    
    rsqplus = 0.045
    rightylabxplus=0.05
    xmargin = 0.005
    cexaxis=1
    cexlab=1
    adj = 0
    
    m <- mask[mask$chr==c & mask$start < poss[2] & mask$end > poss[1],]
    start <- m$start
    end <- m$end

    genes <- dat[dat[,"chr"]==c & dat[,"cdsStart"] < poss[2] & dat[,"cdsEnd"] > poss[1],]
    #par(mar=c(0.5,4,1,4))

    cdsStart <- genes$txStart
    cdsEnd <- genes$txEnd
    geneNames <- genes$name2
    exonCounts <- genes$exonCount

    
    
    plot(c(0,0),c(0,0),type="n",xlim=c(poss[1]-adj,poss[2]+adj),ylim=c(-0.8,0.1),xlab="",yaxt='n',ylab="",main="",xaxt="n", bty="n")

    ord <- order(cdsStart)
    geneNames     <- geneNames[ord]
    
    keep <- !duplicated(geneNames)
    cdsStart    <- cdsStart[ord][keep]
    cdsEnd      <- cdsEnd[ord][keep]
    exonCounts <- exonCounts[ord][keep]
    exonStarts <- genes$exonStarts[ord][keep]
    exonEnds <- genes$exonEnds[ord][keep]
    geneNames     <- geneNames[keep]
    ord <- ord[keep]
    he       <- rep(c(0,-0.18,-0.36,-0.54,-0.72),100)[1:length(geneNames)]-0.05

    if(length(cdsStart)>0){
        segments(cdsStart, he, cdsEnd, he)
        if(sum(keep) < 50) sapply(1:sum(keep),function(x){text((cdsEnd[x]+cdsStart[x])/2,he[x]+0.08,bquote(italic(.(geneNames[x]))),cex=0.7, xpd=NA)})
     
        estart = as.numeric(unlist(sapply(exonStarts,function(y){strsplit(y,",")[[1]]})))
        eend = as.numeric(unlist(sapply(exonEnds,function(y){strsplit(y,",")[[1]]})))
        rect(estart,rep(he,exonCounts)-0.01,eend, rep(he, exonCounts) + 0.01,col="black")
    }
  
    par(mar=c(2.1, 4.2, 0, 2.1))
    fst_col <- "#00A4CCFF"
    fd_col <- "#F95700FF"
    
    plot(x=res$midPos[res$chr==c & (res$start < poss[2] & res$end > poss[1])],y=res$Fst01[res$chr==c & (res$start < poss[2] & res$end > poss[1])], type="l", xlab="", ylab="Window statistic", ylim=c(0,1), xlim=poss, bty="n", col=fst_col,cex.axis=1, cex.lab=1)#col=brewer.pal(7, "Set1")[6])

    rect(xleft=start,xright=end,ybottom=-1.5,ytop=1.5,col=alpha("grey",0.25), border =NA, cex=2)
    abline(h=q_fst, lty=2, col=fst_col)

    
    lines(x=res$midPos[res$chr==c & (res$start < poss[2] & res$end > poss[1])],y=res$fd[res$chr==c & (res$start < poss[2] & res$end > poss[1])], col=fd_col)#col=brewer.pal(7, "Set1")[6])
    abline(h=q_fd, lty=2, col=fd_col)
    mtext(1,text=paste("Position on chromosome ",c,sep=""),line=3,cex=0.7)
    
}



plotSegment("chr5", 72000001, winplot=1e6)
plotSegment("chr7", 22500001, winplot=15e6)

text("MHC locus", x=22.5e6, y=1.63, xpd=NA, cex=0.7)
segments(x0=c(20e6, 25e6), x1=c(20e6, 25e6), y0=1.05, y1=1.6, lty=2, xpd=NA)

text(x=mean(c(26562801, 26700415)), xpd=NA, y=1.125, labels="MILP", cex=0.7)
text(x=mean(c(19289384, 19345780)), xpd=NA, cex=0.7, y=1.49, labels="GPLD1")


fst_col <- "#00A4CCFF"
fd_col <- "#F95700FF"

legend(x=0.9e7, y=-0.6, legend=c(expression("Window F"["st"]), expression("Window f"["d"])), col=c(fst_col, fd_col), lty=1, xpd=NA, bty="n", cex=1)
legend(x=1.5e7, y=-0.6, legend=c(expression("00.5 % quantile F"["st"]), expression("99.9 % quantile f"["d"])), col=c(fst_col, fd_col), lty=2, xpd=NA, bty="n", cex=0.9)

text(x=c(-0.5e6, 2.6e7, -0.7e6), y=c(3.9, 3.9, 1.8),cex=2, labels=c("A", "B", "C"), xpd=NA)
