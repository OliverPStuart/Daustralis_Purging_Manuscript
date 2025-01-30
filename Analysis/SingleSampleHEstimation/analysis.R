
### Script to analyse and plot H from lcWGS

# Libraries

library(dplyr)
library(magrittr)
library(ggplot2)
library(ggh4x)

# Set environment

source("../../config.R")
source(paste0(DATA_DIR,"/colours.R"))
setwd(paste0(WORKING_DIR,"/SingleSampleHEstimation"))

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

# Plot the H versus depth for the two depth classes

# Calculate mean genome-wide H across scaffolds

mean_het <- het %>% group_by(Sample,Depth) %>% 
  mutate(weight_het = Het * Sites) %>%
  dplyr::summarise(mean_het = sum(weight_het)/sum(Sites),
                   Pop=first(Pop),Sex=first(Sex),
                   mean_cov=mean(Coverage),
                   Sites=sum(Sites))
mean_het$Depth <- recode(mean_het$Depth,low = "2-4X",mid = "5-8X")

# Remove wild, as these are outliers for both depth and heterozygosity
# This will create spurious correlation
# Also remove C01220, since her heterozygosity is produced by a different process

p <- mean_het %>% 
  filter(Pop != "Wild",Sample !="C01220") %>%
  ggplot(aes(x=Sites,y=mean_het,fill=Pop,shape=Sex)) + 
  geom_smooth(method="lm",inherit.aes=F,aes(x=Sites,y=mean_het),colour="black") + 
  geom_point(size=3) + 
  scale_shape_manual(values=c(21,24)) +
  my_fill_3 +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill="lightgrey"))) +
  facet_wrap(~Depth,scales="free_x") + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.spacing.x = unit(0.5,"cm")) + 
  labs(x="Sites sampled",y="Genomic H") + 
  scale_y_continuous(labels=c(expression(3%*%10^-5),expression(5%*%10^-5),expression(7%*%10^-5),expression(9%*%10^-5)),
                     breaks=c(3e-5,5e-5,7e-5,9e-5)) + 
  facetted_pos_scales(
    x = list(
      Depth == "2-4X" ~ scale_x_continuous(breaks = c(5e8,5.5e8,6e8,6.5e8),
                                           labels = c(expression(5%*%10^8),expression(5.5%*%10^8),
                                                      expression(6%*%10^8),expression(6.5%*%10^8))),
      Depth == "5-8X" ~ scale_x_continuous(breaks = c(2e8,3e8,4e8,5e8),
                                           labels = c(expression(2%*%10^8),expression(3%*%10^8),
                                                      expression(4%*%10^8),expression(5%*%10^8)))
    )
  ) ; p

# Save figure

png(paste0(FIGURE_DIR,"/het_v_depth_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=7,height=6,units='in')
plot(p)
dev.off()

# Run a simple model to estimate effect of mean_cov on mean_het

mean_het %>% 
  filter(
    Pop != "Wild",
    Depth == "2-4X",
    Sample != "C01220"
  ) %>%
  lm(mean_het ~ Sites,data=.) %>% summary
mean_het %>% 
  filter(
    Pop != "Wild",
    Depth == "5-8X",
    Sample != "C01220"
  ) %>%
  lm(mean_het ~ Sites,data=.) %>% summary

# Refer to model outputs in figure caption

# Calculate mean genome-wide H across scaffolds

mean_het <- het %>% group_by(Sample,Depth) %>% 
  mutate(weight_het = Het * Sites) %>%
  dplyr::summarise(mean_het = sum(weight_het)/sum(Sites),
                   Pop=first(Pop),Sex=first(Sex),
                   Sites=sum(Sites))

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

# Also save table for paper

mean_het %>%
  arrange(Depth,Pop,Sample) %>%
  mutate(Depth=recode(Depth,"low"="Low","mid"="Mid"),
         Pop=recode(Pop,"Wild"="Wild","LHIP"="Hybrid","LHISI"="Captive")) %>%
  select(Sample,Depth,Pop,Sites,mean_het,Sex) %>%
  write.table(.,"mean_het_summarised.txt",col.names=T,row.names=F,sep="\t",quote=F)
