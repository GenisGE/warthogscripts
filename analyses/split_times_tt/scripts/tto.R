# Rscriot tto.R <path to unfolded 2dsfs> <prefix for output files> <assumed mutation rate> <assumed generation time>

args <- commandArgs(trailingOnly=T)

insfs <- args[1]
outprefix <- args[2]
mu <- as.numeric(args[3])
g <- as.numeric(args[4])

dat <- read.table(insfs,stringsAsFactors=F)

sfs <- as.numeric(dat$V2)

# this is just to reformat sfs bin labels so eg 0/0_0/0 -> 00; 0/1_0/0 -> 10, 0/1_1/1 -> 12... 
names(sfs) <- sapply(lapply(strsplit(dat$V1,"_"), function(x) lapply(strsplit(x, split="/"), function(x) sum(as.integer(x)))), paste, collapse="")

# test of treenness to check that outgorup is a proepr oturput. should be 0
y1 <- (2*sfs["10"] + sfs["11"]) / (2*sfs["01"] + sfs["11"]) - (2*sfs["12"] + sfs["11"]) / (2*sfs["21"] + sfs["11"])
y2 <- ((sfs["10"] - sfs["01"]) + 2 * (sfs["20"] - sfs["02"]) + (sfs["21"] - sfs["12"])) / sum(sfs)

outfile <- paste0(outprefix,".treenness.test")
cat(file=outfile,
    "Test of treeness (should be 0 if your close outgroup is a real outgroup):\n\nY1: ",y1,"\nY2:" ,y2, "\n")


# alphas will be used later
alpha1 <- 2 * (sfs["10"] + sfs["12"] + sfs["11"]) / (2*(sfs["10"] + sfs["20"] + sfs["21"]) + sfs["11"])
alpha2 <- 2 * (sfs["01"] + sfs["21"] + sfs["11"]) / (2*(sfs["01"] + sfs["02"] + sfs["12"]) + sfs["11"])


t2a <- (3/2) * ((2*sfs["21"] + sfs["11"]) / alpha2 - sfs["11"]/(alpha1*alpha2)) / sum(sfs)
t2b <- (3/2) * ((2*sfs["12"] + sfs["11"]) / alpha2 - sfs["11"]/(alpha1*alpha2)) / sum(sfs)
t2 <- mean(c(t2a,t2b))

t3a <- (5/2 * (sfs["11"] / (alpha1*alpha2)) - ((2*sfs["21"] + sfs["11"]) / alpha2)) / sum(sfs)
t3b <- (5/2 * (sfs["11"] / (alpha1*alpha2)) - ((2*sfs["12"] + sfs["11"]) / alpha2)) / sum(sfs)
t3 <- mean(c(t3a, t3b))

b1 <- (sfs["10"]/2 + sfs["20"] + sfs["21"] / 2 - ((5 - alpha1 * alpha2) / (alpha1 * alpha2 )) * (sfs["11"] / 4)) / sum(sfs)
b2 <- (sfs["01"]/2 + sfs["02"] + sfs["12"] / 2 - ((5 - alpha1 * alpha2) / (alpha1 * alpha2 )) * (sfs["11"] / 4)) / sum(sfs)

v1 <- NA # TO DO (when we figure out if we care)
v2 <- NA # TO DO (when we figure out if we care)

t4 <- (3/2) * ((t3^2)/t2)

T1 <- b1 - t4
T2 <- b2 - t4

res <- c("alpha1" = unname(alpha1), "alpha2" = unname(alpha2),
         "t2a" = unname(t2a), "t2b" = unname(t2b), "t2" = unname(t2),
         "t3a" = unname(t3a), "t3b" = unname(t3b), "t3" = unname(t3),
         "t4" = unname(t4),
         "v1" = unname(v2), "v2" = unname(v2),
         "b1" = unname(b1), "b2" = unname(b2),
         "T1" = unname(T1), "T2" = unname(T2)
         )

outfile <- paste0(outprefix, ".tto.params.res")
write.table(x = t(t(res)), file = outfile,col.names=F, row.names=T,quote=F)

# scale time in years
res_scaled<- c(
    T1 = unname(res["T1"] * g / mu),
    T2 = unname(res["T2"] * g / mu),
    `mu (assumed)` = mu,
    `g (assumed)` = g
    
)


outfile <- paste0(outprefix, ".tto.split.years")
write.table(x = t(t(res_scaled)), file = outfile,col.names=F, row.names=T,quote=F)
