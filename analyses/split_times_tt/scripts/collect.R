

indir <- "/home/genis/warthogs/analyses/split_times_tt/results/estimates/tt"
outf <- "tt_warthogs.res"


f <- list.files(indir, pattern=".split.years", full.names=T)


res <- lapply(f,read.table,nrows=2)
names(res) <- gsub(".split.years", "", basename(f))

df <- do.call(rbind,res)


outpng <- "/home/genis/warthogs/analyses/split_times_tt/plots/warts_tt1.png"
bitmap(outpng, res=300,w=8,h=6)
par(oma = c(6,0,0,0))
barplot(df$V2/1000000, las=2,
        names.arg=paste(gsub(".[12]$", "", rownames(df)), df$V1, sep="_"),
        ylab="Mya", main="Warthog split times point estimates by TT method")
dev.off()
