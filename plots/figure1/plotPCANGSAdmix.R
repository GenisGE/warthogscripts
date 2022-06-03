source("/home/genis/warthogs/info/loadPopInfo.R")
source("/home/genis/software/evalAdmix/visFuns.R")

fcov <- "/home/jfl323/paper/pca/35CommonNoRepDepthMapHWE_filtered.cov"
fq <- "/home/genis/warthogs/analyses/ngsadmix/results/5/admixResultCommonWart.5.9.qopt_conv"

popf <- "/home/jfl323/paper/filelists/35CommonFromThesis_withPop.filelist"

pop <- read.table(popf,h=F, stringsAsFactors=F)$V1


cov <- read.table(fcov)
ev <- eigen(cov)
e <- ev$vectors
vars <- ev$values/sum(ev$values) * 100



refpops <- c("Namibia", "Ghana", "Tanzania", "Zimbabwe", "Kenya")
popord <- popord[! popord == "Desert"]

q <- as.matrix(read.table(fq))

outpng <- "/home/genis/warthogs/paper_plots/pop_structure/commonwarts_popstructure.png"
bitmap(outpng, h=10, w=5, res=300)
par(mfrow=c(2,1), oma=c(0,1,0,0), xpd=NA)
plot(e[,1], e[,2], pch=21, cex=2, bg=wartCols[pop], cex.lab=1.5,
     xlab=paste0("PC 1 (",round(vars[1], 2),"%)"), ylab=paste0("PC 2 (",round(vars[2], 2),"%)"),cex.lab=1.6)
legend(x=-0.45, y=0.2,
       legend=popord, pt.bg=wartCols[popord],
       pch=21,bty='n', cex=2)


par(mar=c(5,4,4,4) +0.1, xpd=F)
ord <- orderInds(pop=pop, popord=popord)
plotAdmix(q[, orderK(q, refpops=refpops, pop=pop)], ord=ord,pop=pop, colorpal=wartCols[refpops],
          cex.labx=1.5, cex.laby=1.5, rotatelab=30,padj=0.02, main="")
text(x=39, y=0.5, labels="K = 5", cex = 2, xpd=NA)
dev.off()





# TO MAKE PLOT WIHT DESERT WARTHOG
if(FALSE){
fcov <- "/home/jfl323/paper/pca/39inds/39SamplesNoRepDepthMapHWE_filtered.cov"

pop <- read.table("/home/jfl323/paper/filelists/allFromThesis_withPop.filelist",h=F, stringsAsFactors=F)$V1

cov <- read.table(fcov)
ev <- eigen(cov)
e <- ev$vectors
vars <- ev$values/sum(ev$values) * 100



bitmap("warts_popstructure.png", h=10, w=5, res=300)
par(mfrow=c(2,1))
plot(e[,1], e[,2], pch=21, cex=2, bg=wartCols[pop], cex.lab=1.5,
     xlab=paste0("PC 1 (",round(vars[1], 2),"%)"), ylab=paste0("PC 2 (",round(vars[2], 2),"%)"),cex.lab=1.6)
legend(x=-0.5, y=-0.2,
       legend=names(wartCols), pt.bg=wartCols,
       pch=21,bty='n', cex=2)


refpops <- c("Desert", "Namibia", "Ghana", "Tanzania", "Zimbabwe", "Kenya")

fq <- "/home/jfl323/paper/admix/ngsadmix/converged/admix_result.6.13.qopt_conv"

par(mar=c(5,4,4,4) +0.1)
q <- as.matrix(read.table(fq))
ord <- orderInds(pop=pop, popord=popord)
plotAdmix(q[, orderK(q, refpops=refpops, pop=pop)], ord=ord,pop=pop, colorpal=wartCols[refpops],
          cex.labx=1.5, cex.laby=1.5, rotatelab=30,padj=0.02, main="")
text(x=43, y=0.5, labels="K = 6", cex = 2, xpd=NA)

dev.off()




resfold <- "/home/jfl323/paper/admix/ngsadmix/converged"
# use later to make supplementary multiadmixture plot
par(mfrow=c(4,1))
for(k in 3:6){
    fq <- list.files(resfold,pattern=paste0(k,".[1-9]*.qopt"), full.names=T)
    

    q <- as.matrix(read.table(fq))

    ord <- orderInds(pop=pop, popord=popord)

  #  plotAdmix(q[, orderK(q, refpops=refpops, pop=pop)], ord=ord,pop=pop, colorpal=wartCols[refpops])
}


}
