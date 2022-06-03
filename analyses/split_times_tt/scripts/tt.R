# Rscriot tt.R <path to unfolded 2dsfs> <prefix for output files> <assumed mutation rate> <assumed generation time>

args <- commandArgs(trailingOnly=T)

insfs <- args[1]
outprefix <- args[2]
mu <- as.numeric(args[3])
g <- as.numeric(args[4])

dat <- read.table(insfs,stringsAsFactors=F)

sfs <- dat$V2

# this is just to reformat sfs bin labels so eg 0/0_0/0 -> 00; 0/1_0/0 -> 10, 0/1_1/1 -> 12... 
names(sfs) <- sapply(lapply(strsplit(dat$V1,"_"), function(x) lapply(strsplit(x, split="/"), function(x) sum(as.integer(x)))), paste, collapse="")

# sfs needs to be 9 length named vector with count of each bin as entries and names are 00 01 10 02 20 11 12 21 22 

# apply formulas in equation 3 from sjodin et al 2021
# t1 and t2 are divergence time in generations for pop 1 and pop2 from ancstral. Ideally they should be the same assuming smae mutation rate in each branch and no messy sample
# the rest of the parameters I am not sure what they are I just copy the formula, maybe should read the paper better
res <- c(
    alpha1 = unname((2 * sfs["11"]) / (2 * sfs["21"] + sfs["11"])),
    alpha2 = unname((2 * sfs["11"]) / (2 * sfs["12"] + sfs["11"])),
    theta = unname(3/8 * (((2*sfs["21"] + sfs["11"]) * (2 * sfs["12"] + sfs["11"])) / (sfs["11"])) / sum(sfs)),
    t1 = unname((sfs["10"]/2 + sfs["20"] - (((2*sfs["21"] + sfs["11"]) * (6*sfs["12"] + sfs["11"])) / (8*sfs["11"]))) / sum(sfs)),
    t2 = unname((sfs["01"]/2 + sfs["02"] - (((6*sfs["21"] + sfs["11"]) * (2*sfs["12"] + sfs["11"])) / (8*sfs["11"]))) / sum(sfs)),
    v1 = unname(((sfs["10"] + sfs["12"]) / 2 - sfs["11"] * ((2*sfs["20"] - sfs["12"]) / (2*sfs["21"] - sfs["11"]))) / sum(sfs)),
    v2 = unname(((sfs["01"] + sfs["21"]) / 2 - sfs["11"] * ((2*sfs["02"] - sfs["21"]) / (2*sfs["12"] - sfs["11"]))) / sum(sfs))
)


outfile <- paste0(outprefix, ".tt.params.res")
write.table(x = t(t(res)), file = outfile,col.names=F, row.names=T,quote=F)


# scale time in years
res_scaled<- c(
    T1 = unname(res["t1"] * g / mu),
    T2 = unname(res["t2"] * g / mu),
    Na = unname(res["theta"] / mu),
    `mu (assumed)` = mu,
    `g (assumed)` = g
)

outfile <- paste0(outprefix, ".tt.split.years")
write.table(x = t(t(res_scaled)), file = outfile ,col.names=F, row.names=T,quote=F)
