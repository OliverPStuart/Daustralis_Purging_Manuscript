### Script to analyse and plot lcWGS PCA

# Set environment

source("../../config.R")
setwd(paste0(WORKING_DIR,"/PCangsd"))

# Libraries

library(dplyr)
library(magrittr)
library(ggplot2)
library(stringi)

# My favourite function

`%ni%` <- Negate(`%in%`)

# Get colours

source(paste0(DATA_DIR,"/colours.R"))

# Get samples

samples <- data.frame(Sample=stri_reverse(gsub(x=unlist(lapply(strsplit(stri_reverse(readLines("all_bam_list.txt")),
                                                                        split="/"),`[[`,1)),pattern="mab\\.","")),
                      Index=1:19)


details <- read.table(paste0(REF_DIR,"/wgs_sample_details.txt"),
                      stringsAsFactors=F,sep="\t",header=T)
samples <- merge(samples,details)
samples <- samples[order(samples$Index),]

# Rename populations

samples$Pop <- as.factor(samples$Pop)

# Load file

file = "pcangsd_all_thin_maf_10.cov"
cov <- as.matrix(read.table(file))

# Calculate eigenvectors and add to samples data.frame

PCs <- eigen(cov)
samples$PC1 <- PCs$vectors[,1]
samples$PC2 <- PCs$vectors[,2]

# Extract components loadings

loadings <- round(PCs$values[1:2] / sum(PCs$values),2) * 100

# Plot

p <- ggplot(samples,aes(x=PC1,y=PC2,fill=Pop)) + 
  geom_point(size=4,pch=21) + 
  my_fill_3 + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.grid = element_blank()) + 
  labs(x=paste0("PC1 (",loadings[1],"%)"),
       y=paste0("PC2 (",loadings[2],"%)"))

# Save figure

png(paste0(FIGURE_DIR,"/pcangsd_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=5,height=5,units='in')
plot(p)
dev.off()

# Now do the same but for all of the maf filters
# Does it make a difference?

file = "pcangsd_all_thin_maf.cov"
cov <- as.matrix(read.table(file))
PCs <- eigen(cov)
samples$PC1 <- PCs$vectors[,1]
samples$PC2 <- PCs$vectors[,2]

pc_data <- data.frame(Sample = samples$Sample,
                      Index = samples$Index,
                      Sex = samples$Sex,
                      Pop = samples$Pop,
                      PC1 = PCs$vectors[,1],
                      PC2 = PCs$vectors[,2],
                      MAF = "0.00")

file = "pcangsd_all_thin_maf_05.cov"
cov <- as.matrix(read.table(file))
PCs <- eigen(cov)
samples$PC1 <- PCs$vectors[,1]
samples$PC2 <- PCs$vectors[,2]

pc_data <- rbind(pc_data,
                 data.frame(Sample = samples$Sample,
                      Index = samples$Index,
                      Sex = samples$Sex,
                      Pop = samples$Pop,
                      PC1 = PCs$vectors[,1],
                      PC2 = PCs$vectors[,2],
                      MAF = "0.05"))

file = "pcangsd_all_thin_maf_10.cov"
cov <- as.matrix(read.table(file))
PCs <- eigen(cov)
samples$PC1 <- PCs$vectors[,1]
samples$PC2 <- PCs$vectors[,2]

pc_data <- rbind(pc_data,
                 data.frame(Sample = samples$Sample,
                            Index = samples$Index,
                            Sex = samples$Sex,
                            Pop = samples$Pop,
                            PC1 = PCs$vectors[,1],
                            PC2 = PCs$vectors[,2],
                            MAF = "0.10"))

file = "pcangsd_all_thin_maf_15.cov"
cov <- as.matrix(read.table(file))
PCs <- eigen(cov)
samples$PC1 <- PCs$vectors[,1]
samples$PC2 <- PCs$vectors[,2]

pc_data <- rbind(pc_data,
                 data.frame(Sample = samples$Sample,
                            Index = samples$Index,
                            Sex = samples$Sex,
                            Pop = samples$Pop,
                            PC1 = PCs$vectors[,1],
                            PC2 = PCs$vectors[,2],
                            MAF = "0.15"))

p <- ggplot(pc_data,aes(x=PC1,y=PC2,fill=Pop)) + 
  geom_point(size=4,pch=21) + 
  my_fill_3 + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.grid = element_blank()) + 
  facet_wrap(~MAF)

# Save figure

png(paste0(FIGURE_DIR,"/pcangsd_maf_filter_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=7,height=7,units='in')
plot(p)
dev.off()



