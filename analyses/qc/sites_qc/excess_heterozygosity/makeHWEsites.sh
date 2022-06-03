python ~/software/pcangsd/pcangsd.py -beagle ../beagle/filteredBeagleForHWE/35ComFromThesisAutosomes.beagle.gz -e 3 -inbreedSites -sites_save -inbreed 1 -o genHWEsites
python ~/software/pcangsd/pcangsd.py -beagle ../beagle/filteredBeagleForHWE/35ComFromThesisAutosomes.beagle.gz -e 3 -hwe genHWEsites.lrt.sites.npy -sites_save -inbreed 1 -o HWEsitesRemoved

