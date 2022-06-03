

args <- commandArgs(trailingOnly=T)

indir <- args[1]
outpng <- args[2]



pattern <- ".tt.split.years"
plottitle <- "Ancestral population size point estimates based on TT method"


f <- list.files(indir, pattern=pattern, full.names=T)


res <- lapply(f,read.table,nrows=1, skip=2)
names(res) <- gsub(pattern, "", basename(f))

df <- do.call(rbind,res)


bitmap(outpng, res=300, w=8, h=6)
par(oma = c(6,0,0,0))
barplot(df$V2/1000, las=2,
        names.arg = rownames(df),
        ylab="Population size (10^3)", main=plottitle)
dev.off()

