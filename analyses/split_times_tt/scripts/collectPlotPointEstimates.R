
args <- commandArgs(trailingOnly=T)

indir <- args[1]
outpng <- args[2]

method <- "tt"
if(length(args)>2) method <- args[3]

if(method == "tt"){
    pattern <- ".tt.split.years"
    plottitle <- "Split times point estimates based on TT method"
}else if(method == "tto"){
    pattern <- ".tto.split.years"
    plottitle <- "Split times point estimates based on TTo method"
    
}else{
    stop("Method (3rd argument) must be either 'tt' 'tto' or nothing")
}

    
f <- list.files(indir, pattern=pattern, full.names=T)


res <- lapply(f,read.table,nrows=2)
names(res) <- gsub(pattern, "", basename(f))

df <- do.call(rbind,res)


bitmap(outpng, res=300, w=8, h=6)
par(oma = c(6,0,0,0))
barplot(df$V2/1000000, las=2,
        names.arg=paste(gsub(".[12]$", "", rownames(df)), df$V1, sep="_"),
        ylab="Mya", main=plottitle)
dev.off()

