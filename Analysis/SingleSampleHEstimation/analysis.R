
### Script to analyse and plot H from lcWGS

# Set environment

source("../../config.R")
source(paste0(DATA_DIR,"/colours.R"))
setwd(paste0(WORKING_DIR,"/SingleSampleHEstimation"))

# Libraries

library(dplyr)
library(magrittr)
library(ggplot2)

# My favourite function

`%ni%` <- Negate(`%in%`)

# Load data

# Get sample information, heterozygosity estimates, and depth

het <- read.table("het_estimates.txt",
                  stringsAsFactors=F,sep="\t",header=T)
samples <- read.table(paste0(REF_DIR,"/wgs_sample_details.txt"),
                      stringsAsFactors=F,sep="\t",header=T)

# In case you're using the old snakefile which doesn't exclude the low coverage
# This line will remove them

lc <- c("C01216","C01233","C01218","C01226","C10133")
samples <- samples[samples$Sample %ni% lc,]

# Merge everything into one table 

het <- merge(het,samples,by="Sample")

# Also get the lengths of scaffolds for some QC

lengths <- read.table(paste0(REF_DIR,"/autosome_regions.bed"),
                      header=F,stringsAsFactors=F,sep="\t")[,c(1,3)]
lengths$V3 <- lengths$V3+1
colnames(lengths) <- c("Chrom","Length")
het <- merge(het,lengths,by="Chrom")

# Now calculate the proportion of sites covered by each sample for each scaffold
# i.e. What is our sample size for estimating heterozygosity?

het$Coverage <- het$Sites/het$Length

# Calculate mean genome-wide H across scaffolds

mean_het <- het %>% group_by(Sample,Depth) %>% 
  mutate(weight_het = Het * Sites) %>%
  dplyr::summarise(mean_het = sum(weight_het)/sum(Sites),
                   Pop=first(Pop),Sex=first(Sex))

# Rename populations

mean_het$Pop <- as.factor(mean_het$Pop)
mean_het$Pop <- ordered(mean_het$Pop,
                        levels = c("Wild", "LHIP", "LHISI"))

# Now plot

p <- mean_het %>% filter(Depth == "mid") %>% 
  ggplot() + 
  stat_summary(aes(x=Pop,y=mean_het, group=Pop),
               fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5) +
  geom_point(aes(x=Pop,y=mean_het, fill=Pop),pch=21,colour="black",
             size=4,position=position_jitterdodge(dodge.width=0.5,jitter.width=0.5)) +
  my_fill_3 + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.text.x=element_blank(),axis.title.x=element_blank(),
        panel.grid=element_blank()) +
  ylab("Genomic H") +
  scale_y_continuous(breaks=c(0,5e-5,1e-4,1.5e-4),
                     labels = c(expression(0),
                                expression(0.5%*%10^-4),
                                expression(1%*%10^-4),
                                expression(1.5%*%10^-4)))

# Save figure

png(paste0(FIGURE_DIR,"/mean_het_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=5,height=5,units='in')
plot(p)
dev.off()

# Test difference
# We test low and mid, and excluding/including C01220 outlier

mean_het %>% filter(Pop != "Wild",Depth == "low") %>% t.test(mean_het~Pop,data=.)
mean_het %>% filter(Pop != "Wild",Depth == "low",Sample != "C01220") %>% t.test(mean_het~Pop,data=.)

mean_het %>% filter(Pop != "Wild",Depth == "mid") %>% t.test(mean_het~Pop,data=.)
mean_het %>% filter(Pop != "Wild",Depth == "mid",Sample != "C01220") %>% t.test(mean_het~Pop,data=.)
