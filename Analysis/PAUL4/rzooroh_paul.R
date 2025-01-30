
setwd("/Volumes/Alter/Daus_WGS_Paper/Analysis/PAUL4/")

file="paul4_gt_rzooroh.txt"

library(RZooRoH)

gt <- zoodata(genofile = file, min_maf=0,poscol=2,chrcol=1,supcol=4)
