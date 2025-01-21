### Script to analyse and plot H from lcWGS

setwd("/Volumes/Alter/Daus_WGS_Paper/Analysis/DeleteriousMutations/")

############################################
###           Environment setup          ###
############################################

# Libraries

library(ggplot2)
library(dplyr)
library(ggnewscale)
library(ggResidpanel)
library(reshape)
library(Hmisc)
library(MASS)
library(foreign)
library(DescTools)
library(brms)
library(ggridges)

# Set environment

source("../../config.R")
setwd(paste0(WORKING_DIR,"/DeleteriousMutations"))

# Functions

`%ni%` <- Negate(`%in%`)
fold <- function(x) {
  ifelse(x > 0.5, 1 - x, x)
}

# Get colours

source(paste0(DATA_DIR,"/colours.R"))

############################################
###            Data preparation          ###
############################################

lhisi <- read.table("../AlleleFrequencies/lhisi.mafs.gz",
                    header=T,stringsAsFactors=F,sep="\t")[,c(1,2,6,7)]
lhip <- read.table("../AlleleFrequencies/lhip.mafs.gz",
                   header=T,stringsAsFactors=F,sep="\t")[,c(1,2,6,7)]
wild <- read.table("../AlleleFrequencies/wild.mafs.gz",header=T,
                   stringsAsFactors=F,sep="\t")[,c(1,2,6,7)]

colnames(lhisi) <- c("Scaffold","Position","MAF_lhisi","N_lhisi")
colnames(lhip) <- c("Scaffold","Position","MAF_lhip","N_lhip")
colnames(wild) <- c("Scaffold","Position","MAF_wild","N_wild")

# Read in effect data

effects <- read.table("vars_classified.txt",
                      header=T,stringsAsFactors=F,sep="\t")[,c(1,2,3,4,5,6,7,8,10)]

# Get site information

scafs <- readLines(paste0(REF_DIR,"/autosome_names"))

sites <- data.frame(Scaffold=character(),
                    Position=character(),
                    Ref=character(),
                    Alt=character())

for(i in 1:length(scafs)){
  tmp <- eval(parse(text=paste0("read.table('../AlleleFrequencies/",
                                scafs[i],
                                "_sites',header=F,stringsAsFactors=F)")))
  sites <- rbind(sites,tmp)
}

rm(tmp) ; colnames(sites) <- c("Scaffold","Position","Ref","Alt")

# Give all sites their alleles

effects <- merge(effects,sites)
rm(sites)

# Give variants inside genes their AED score

scores <- read.table("gene_scores.txt",header=T,stringsAsFactors = F)

intergenic <- effects %>% filter(is.na(Gene) == T)
genic <- effects %>% filter(is.na(Gene) == F)
genic <- merge(genic,scores)
effects <- plyr::rbind.fill(genic,intergenic)
rm(genic,intergenic)

# For each population, score every site by that population's frequency for the allele

classify_maf <- function(x){
  if(x <= 0.05){return("absent")} else
    if (x & x < 0.95 ){return("segregating")} else
      if (x >= 0.95 ){return("fixed")}
}

wild$fate_wild <- sapply(wild$MAF_wild,classify_maf)
lhisi$fate_lhisi <- sapply(lhisi$MAF_lhisi,classify_maf)
lhip$fate_lhip <- sapply(lhip$MAF_lhip,classify_maf)

# Now get all of these into the effects dataset

effects <- merge(effects,wild,all.x=T)
effects <- merge(effects,lhisi,all.x=T)
effects <- merge(effects,lhip,all.x=T)

# Make an ID column for later

effects$ID <- paste0(effects$Scaffold,"_",effects$Position)

rm(lhisi,lhip,wild)

# Tally up numbers of segregating variants in all groups

# How many sites detected?
effects %>% filter(!is.na(fate_wild)) %>% nrow()
effects %>% filter(!is.na(fate_lhisi)) %>% nrow()
effects %>% filter(!is.na(fate_lhip)) %>% nrow()

# How many segregating of each kind?
WILD <- effects %>% filter(fate_wild=="segregating") %>% xtabs(data=.,formula=~Variant_Effect)
LHISI <- effects %>% filter(fate_lhisi=="segregating") %>% xtabs(data=.,formula=~Variant_Effect)
LHIP <- effects %>% filter(fate_lhip=="segregating") %>% xtabs(data=.,formula=~Variant_Effect)
WILD
LHISI
LHIP

# Chi-squared test to see if the proportion are different
rbind(WILD,LHISI,LHIP) %>%
  chisq.test()
WILD/sum(WILD)
LHIP/sum(LHIP)
LHISI/sum(LHISI)

# Classifying variants based on presence in putative ROH in the three populations

lhip <- read.table("LHIP.cov",header=F,stringsAsFactors=F,sep="\t")
lhisi <- read.table("LHISI.cov",header=F,stringsAsFactors=F,sep="\t")
wild <- read.table("Wild.cov",header=F,stringsAsFactors=F,sep="\t")
colnames(lhip) <- c("Scaffold","Position","Coverage_lhip")
colnames(lhisi) <- c("Scaffold","Position","Coverage_lhisi")
colnames(wild) <- c("Scaffold","Position","Coverage_wild")

tmp <- merge(lhip,lhisi)
tmp <- merge(tmp,wild)
 
effects <- merge(effects,tmp,all.x=T)
 
rm(tmp,lhisi,lhip,wild)

# Some refactorisation to order specific variables

effects$Variant_Effect <- factor(effects$Variant_Effect,
                                 levels=c("MODIFIER","LOW","MODERATE","HIGH"),
                                 labels=c("MODIFIER","LOW","MODERATE","HIGH"))

############################################
###             Colour scheme            ###
############################################

new_colours <- read.delim(paste0(REF_DIR,"/two_pop_2d_colour_scale.txt"),
                          header=F,stringsAsFactors = F)
colnames(new_colours) <- c("Effect","Captive","Hybrid","Wild")

new_colours_n <- new_colours %>% tidyr::gather(-Effect,key="Population",value="Colour")

legend_plot <- new_colours_n %>%
  ggplot() + 
  geom_point(pch=21,size=7,
             aes(x=Population,y=Effect,fill=Effect),
             filter(new_colours_n,Population == "Hybrid")) +
  scale_fill_manual(values=rev(new_colours[,2])) + 
  new_scale_fill() +
  geom_point(pch=21,size=7,
             aes(x=Population,y=Effect,fill=Effect),
             filter(new_colours_n,Population == "Captive")) +
  scale_fill_manual(values=rev(new_colours[,3])) + 
  theme_bw() +
  theme(axis.ticks=element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),legend.position="none",
        axis.title=element_blank(),
        axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=1.1),
        axis.text.y=element_text(size=12)) +
  #  scale_x_discrete(expand=c(-0.1,-0.1)) + 
  scale_y_discrete(
    #expand=c(-0.1,-0.1),
    labels=c("High","Moderate","Low","Non-coding"))

# Plot

png(paste0(FIGURE_DIR,"/frequency_change_legend_",format(Sys.time(),"%Y%m%d"),".png"),
    width=1.7,height=2,units='in',res=300)
plot(legend_plot)
dev.off()

############################################
###           Raw AF difference         ###
############################################

#########
######### NB: No longer used, see Rxy code below
#########

# # Get wild alleles
# 
# wild <- effects %>% filter(fate_wild != "absent" & !is.na(fate_wild))
# 
# # Calculate difference in frequency
# 
# wild$lhisi_change <- wild$MAF_lhisi - wild$MAF_wild
# wild$lhip_change <- wild$MAF_lhip - wild$MAF_wild
# wild <- wild %>% tidyr::drop_na(lhisi_change,lhip_change)
# 
# # Plotting mean+sd difference per variant effect class
# # Excluding those variants in each population which are absent
# 
# tmp_lhisi <- wild %>% group_by(Variant_Effect) %>% 
# #  filter(is.na(AED) | AED > 0.25) %>%
#   filter(fate_lhisi != "absent") %>%
#   dplyr::summarise(mean_lhisi=mean(lhisi_change),
#                    sd_lhisi=sd(lhisi_change)) %>%
#   tidyr::gather(key, value, -Variant_Effect) %>%
#   tidyr::extract(key, c("stat", "pop"), "(.+)_(.+)") %>%
#   tidyr::spread(stat, value)
# tmp_lhip <- wild %>% group_by(Variant_Effect) %>% 
#   filter(fate_lhip != "absent") %>%
# #  filter(is.na(AED) | AED > 0.25) %>%
#   dplyr::summarise(mean_lhip=mean(lhip_change),
#                    sd_lhip=sd(lhip_change)) %>%
#   tidyr::gather(key, value, -Variant_Effect) %>%
#   tidyr::extract(key, c("stat", "pop"), "(.+)_(.+)") %>%
#   tidyr::spread(stat, value)
# 
# change_data <- rbind(tmp_lhisi,tmp_lhip)
# 
# p <- change_data %>% 
#   ggplot() + 
#   scale_x_continuous(limits=c(-0.7,0.1)) + 
#   geom_errorbar(aes(fill=factor(Variant_Effect,levels=c("HIGH","MODERATE","LOW","MODIFIER")),
#                     x=mean,y=pop,xmin=mean-sd,xmax=mean+sd),
#                 width=0,size=0.7,position=position_dodge(width=0.6),
#                 filter(change_data,pop=="lhisi")) + 
#   geom_point(aes(fill=factor(Variant_Effect,levels=c("HIGH","MODERATE","LOW","MODIFIER")),
#                  x=mean,y=pop),
#              size=4.5,pch=21,position=position_dodge(width=0.6),
#              filter(change_data,pop=="lhisi")) + 
#   scale_fill_manual(values=rev(new_colours[,3])) +
#   new_scale_fill() +
#   geom_errorbar(aes(fill=factor(Variant_Effect,levels=c("HIGH","MODERATE","LOW","MODIFIER")),
#                     x=mean,y=pop,xmin=mean-sd,xmax=mean+sd),
#                 width=0,size=0.7,position=position_dodge(width=0.6),
#                 filter(change_data,pop=="lhip")) + 
#   geom_point(aes(fill=factor(Variant_Effect,levels=c("HIGH","MODERATE","LOW","MODIFIER")),
#                  x=mean,y=pop),
#              size=4.5,pch=21,position=position_dodge(width=0.6),
#              filter(change_data,pop=="lhip")) + 
#   scale_fill_manual(values=rev(new_colours[,2])) +
#   geom_vline(xintercept=0,linetype="dashed",colour="black",size=1) + 
#   theme_bw() + 
#   theme(axis.ticks=element_blank(),
#         panel.background=element_blank(),
#         panel.grid.minor=element_blank(),
#         strip.background = element_blank(),
#         strip.text = element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.x=element_text(size=10),
#         axis.text.y=element_blank(),
#         panel.grid=element_blank(),
#         legend.position="none") + 
#   labs(x="Frequency difference\nin captivity") +
#   scale_y_discrete(labels=c("LHIP","LHISI")) ; p
# 
# # Plot
# 
# png(paste0(FIGURE_DIR,"/frequency_change_",format(Sys.time(),"%Y%m%d"),".png"),
#     width=6,height=3,units='in',res=300)
# plot(p)
# dev.off()
# 
# # Now run the tests
# # We exclude WildPresent-CaptiveAbsent alleles here
# # However, we include WildPresent-HybridAbsent
# 
# tmp_lhisi <- wild %>% 
#   filter(fate_lhisi != "absent") %>% 
# #  filter(is.na(AED) | AED > 0.25) %>%
#   dplyr::select(Variant_Effect,lhisi_change) %>% mutate(Pop="LHISI")
# tmp_lhip <- wild %>%
# #  filter(is.na(AED) | AED > 0.25) %>%
#   filter(fate_lhip != "absent") %>% 
#   dplyr::select(Variant_Effect,lhip_change) %>% mutate(Pop="LHIP")
# colnames(tmp_lhip)[2] <- "change" ; colnames(tmp_lhisi)[2] <- "change"
# tmp <- rbind(tmp_lhisi,tmp_lhip)
# 
# # Make factor ordered
# tmp$Variant_Effect <- ordered(tmp$Variant_Effect, 
#                               levels = c("MODIFIER",
#                                          "LOW",
#                                          "MODERATE",
#                                          "HIGH")) 
# 
# # BRMS model? Why discuss significance at all
# 
# lhisi_brms <- brms::brm(change~Variant_Effect,
#                            data = tmp[tmp$Pop=="LHISI",],
#                            family = gaussian(), 
#                            chains = 4,
#                            iter = 3000, warmup = 1000,cores=4,
#                         silent=2)
# summary(lhisi_brms)
# plot(lhisi_brms)
# hypothesis(lhisi_brms,
#            c("Variant_Effect.L<0",
#              "Variant_Effect.Q<0",
#              "Variant_Effect.C<0"),
#            alpha=0.05)
# 
# lhip_brms <- brms::brm(change~Variant_Effect,
#                         data = tmp[tmp$Pop=="LHIP",],
#                         family = gaussian(), 
#                         chains = 4,
#                         iter = 3000, warmup = 1000,cores=4,
#                        silent=2)
# summary(lhip_brms)
# plot(lhip_brms)
# hypothesis(lhip_brms,
#            c("Variant_Effect.L<0",
#              "Variant_Effect.Q<0",
#              "Variant_Effect.C<0"),
#            alpha=0.05)
# 
# # We'll also do thi
# tmp_lhisi <- wild %>% 
#   filter(fate_lhisi != "absent") %>% 
#   #  filter(is.na(AED) | AED > 0.25) %>%
#   dplyr::select(Variant_Effect,lhisi_change) %>% mutate(Pop="LHISI")
# tmp_lhip <- wild %>%
#   filter(fate_lhip != "absent") %>% 
#   #  filter(is.na(AED) | AED > 0.25) %>%
#   dplyr::select(Variant_Effect,lhip_change) %>% mutate(Pop="LHIP")
# colnames(tmp_lhip)[2] <- "change" ; colnames(tmp_lhisi)[2] <- "change"
# tmp <- rbind(tmp_lhisi,tmp_lhip)
# 
# lhisi_brms <- brms::brm(change~Variant_Effect,
#                         data = tmp[tmp$Pop=="LHISI",],
#                         family = gaussian(), 
#                         chains = 4,
#                         iter = 3000, warmup = 1000,cores=4,
#                         silent=2)
# summary(lhisi_brms)
# plot(lhisi_brms)
# hypothesis(lhisi_brms,c("Variant_EffectLOW<0",
#                         "Variant_EffectMODERATE<0",
#                         "Variant_EffectHIGH<0"),alpha=0.05)
# 
# lhip_brms <- brms::brm(change~Variant_Effect,
#                        data = tmp[tmp$Pop=="LHIP",],
#                        family = gaussian(), 
#                        chains = 4,
#                        iter = 3000, warmup = 1000,cores=4,
#                        silent=2)
# summary(lhip_brms)
# plot(lhip_brms)
# hypothesis(lhip_brms,c("Variant_EffectLOW<0",
#                         "Variant_EffectMODERATE<0",
#                         "Variant_EffectHIGH<0"),alpha=0.05)

############################################
###                 Rxy                  ###
############################################

# We do block jackknifing following Do et al. (2015).

# Testing block size shows that block size makes no difference.

# For this, we care about alleles which are present in the focal population. If
# an allele is absent in a population, this could be because selection removed
# it after it was introduced or because it was never introduced into the captive
# breeding program in the first place. This is conservative, because we may
# be removing some strong signals of selection (true positives) for the sake of
# reducing the false positive rate.

# In order to compare differences in RXY, we use the same blocks for all variant
# types, i.e. we define blocks based on all variants, then resample the whole
# dataset and recalculate RXY. This allows us to test if the difference in RXY
# is greater than 0 between groups using a paired t-test.

# LHISI

tmp <- effects %>% 
  filter(fate_lhisi != "absent") %>% 
  mutate(LXY=MAF_lhisi*(1-MAF_wild),LYX=MAF_wild*(1-MAF_lhisi)) %>%
  filter(!is.na(LXY),!is.na(LYX),LXY!=Inf,LYX!=Inf,LXY!=-Inf,LYX!=-Inf) %>%
  dplyr::select(Scaffold,Position,LXY,LYX,Variant_Effect)

# Order the data.frame by Scaffold and Position
tmp <- tmp[order(tmp$Scaffold, tmp$Position), ]

# Initialise block column
tmp$Block <- NA

# Get scaffolds
scafs <- unique(tmp$Scaffold)

# Total number of blocks desired
# This is what gives us exactly 100 blocks of roughly equal size
total_blocks <- 94

# Total number of entries
total_entries <- nrow(tmp)

# Calculate the rough number of entries per block
rough_block_size <- ceiling(total_entries / total_blocks)

# Loop over scaffolds to define blocks algorithmically
for(i in 1:length(scafs)){
  
  # Get entries
  entries <- tmp[tmp$Scaffold==scafs[i],]
  
  # How many blocks roughly in scaffold?
  blocks_in <- ceiling(nrow(entries)/rough_block_size)
  
  # Make a vector of block assignments
  blocks_assign <- sort(rep(1:blocks_in, length.out=nrow(entries)))
  
  # Add block assignments to data.frame
  # If it's not the first scaffold, add previous values to get correct block ID
  if(i == 1){
    
    tmp$Block[tmp$Scaffold==scafs[i]] <- blocks_assign
    
  } else {
    
    blocks_assign <- blocks_assign + max(tmp$Block[tmp$Scaffold==scafs[i-1]])
    tmp$Block[tmp$Scaffold==scafs[i]] <- blocks_assign
    
  }
}

n_blocks <- max(tmp$Block)

# Jackknife, blocks
for(i in 1:n_blocks){
  
  # Remove block, calculate RXY, assign to temp data.frame
  assign(paste0("jack_rxy_",i),
         tmp %>% filter(Block != i) %>% 
           group_by(Variant_Effect) %>%
           dplyr::summarise(RXY=sum(LXY)/sum(LYX),Jack=i))
  
}

# Collate results 
table_names <- ls(pattern = "^jack_rxy_")
table_list <- mget(table_names)
rxy_lhisi <- do.call(rbind, table_list) %>% mutate(Pop="LHISI")
rm(list=ls(pattern = "^jack_rxy_"))

# LHIP

tmp <- effects %>% 
  filter(fate_lhip != "absent") %>% 
  mutate(LXY=MAF_lhip*(1-MAF_wild),LYX=MAF_wild*(1-MAF_lhip)) %>%
  filter(!is.na(LXY),!is.na(LYX),LXY!=Inf,LYX!=Inf,LXY!=-Inf,LYX!=-Inf) %>%
  dplyr::select(Scaffold,Position,LXY,LYX,Variant_Effect)

# Order the data.frame by Scaffold and Position
tmp <- tmp[order(tmp$Scaffold, tmp$Position), ]

# Initialise block column
tmp$Block <- NA

# Get scaffolds
scafs <- unique(tmp$Scaffold)

# Total number of blocks desired
# Can't get exactly 100 equal sized blocks for LHIP, go with 102
total_blocks <- 93

# Total number of entries
total_entries <- nrow(tmp)

# Calculate the rough number of entries per block
rough_block_size <- ceiling(total_entries / total_blocks)

# Loop over scaffolds to define blocks algorithmically
for(i in 1:length(scafs)){
  
  # Get entries
  entries <- tmp[tmp$Scaffold==scafs[i],]
  
  # How many blocks roughly in scaffold?
  blocks_in <- ceiling(nrow(entries)/rough_block_size)
  
  # Make a vector of block assignments
  blocks_assign <- sort(rep(1:blocks_in, length.out=nrow(entries)))
  
  # Add block assignments to data.frame
  # If it's not the first scaffold, add previous values to get correct block ID
  if(i == 1){
    
    tmp$Block[tmp$Scaffold==scafs[i]] <- blocks_assign
    
  } else {
    
    blocks_assign <- blocks_assign + max(tmp$Block[tmp$Scaffold==scafs[i-1]])
    tmp$Block[tmp$Scaffold==scafs[i]] <- blocks_assign
    
  }
}

n_blocks <- max(tmp$Block)

# Jackknife, blocks
for(i in 1:n_blocks){
  
  # Remove block, calculate RXY, assign to temp data.frame
  assign(paste0("jack_rxy_",i),
         tmp %>% filter(Block != i) %>% 
           group_by(Variant_Effect) %>%
           dplyr::summarise(RXY=sum(LXY)/sum(LYX),Jack=i))
  
}

# Collate results 
table_names <- ls(pattern = "^jack_rxy_")
table_list <- mget(table_names)
rxy_lhip <- do.call(rbind, table_list) %>% mutate(Pop="LHIP")
rm(list=ls(pattern = "^jack_rxy_"))

# And combine pops
jackknifes <- rbind(rxy_lhisi,
                     rxy_lhip)
jackknifes$Variant_Effect <- factor(jackknifes$Variant_Effect,
                                    levels=c("HIGH","MODERATE","LOW","MODIFIER"),
                                    labels=c("HIGH","MODERATE","LOW","MODIFIER"))

# Plot the distributions of Rxy for each
jackknife_se <- function(x){
  sqrt((length(x) - 1) / length(x) * sum((x - mean(x))^2))
}

plotting_results <- jackknifes %>%
  group_by(Variant_Effect,Pop) %>%
  dplyr::summarise(Mean=mean(RXY),
                   SE=jackknife_se(RXY),
                   Upper95=Mean+1.96*SE,
                   Lower95=Mean-1.96*SE) %>%
  ungroup

p <- plotting_results %>% 
  ggplot() +
  geom_errorbar(aes(y=Pop,xmax=Mean+SE,xmin=Mean-SE,group=Variant_Effect),
                size=0.75,position=position_dodge(width=0.6),width=0,
                filter(plotting_results,Pop=="LHISI")) + 
  geom_point(aes(y=Pop,x=Mean,fill=Variant_Effect,group=Variant_Effect),
             position=position_dodge(width=0.6),pch=21,size=4.5,
             filter(plotting_results,Pop=="LHISI"))  +
  scale_fill_manual(values=rev(new_colours[,3])) + 
  new_scale_fill() +
  geom_errorbar(aes(y=Pop,xmax=Mean+SE,xmin=Mean-SE,group=Variant_Effect),
                size=0.75,position=position_dodge(width=0.6),width=0,
                filter(plotting_results,Pop=="LHIP")) + 
  geom_point(aes(y=Pop,x=Mean,fill=Variant_Effect),
             position=position_dodge(width=0.6),pch=21,size=4.5,
             filter(plotting_results,Pop=="LHIP")) + 
  scale_fill_manual(values=rev(new_colours[,2])) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        legend.position="none"
        ) + 
  scale_x_continuous(
    limits=c(0.45,0.93),
    breaks=seq(0.5,0.9,0.1)) +
  labs(x=expression(R[XY] ~ "(relative to wild individuals)")) ; p

# Plot

png(paste0(FIGURE_DIR,"/rxy_",format(Sys.time(),"%Y%m%d"),".png"),
    width=6,height=3,units='in',res=300)
plot(p)
dev.off()

# Now we compute the pairwise differences between mutation types within each
# jackknife resample and calculate the t-statistic on the distribution of
# differences. We then compute the probability of observing this t-value under
# the null hypothesis that the RXY values are the same in both types. It is a
# one sided test since we expect the value of RXY to go down to go down each
# time we go up a type (e.g. from low to moderate impact).

# Generate a data.frame of combinations manually then loop over combinations.

combs <- data.frame(Type1=c("MODIFIER","MODIFIER","MODIFIER",
                            "LOW","LOW",
                            "MODERATE"),
                    Type2=c("LOW","MODERATE","HIGH",
                            "MODERATE","HIGH",
                            "HIGH")) %>% arrange(Type2)
pops <- c("LHISI","LHIP")

for(pop in pops){
  for(i in 1:nrow(combs)){
    
    # Get vectors of values
    Type1=jackknifes %>% arrange(Jack) %>% filter(Variant_Effect == combs$Type1[i],Pop==pop) %>% pull(RXY)
    Type2=jackknifes %>% arrange(Jack) %>% filter(Variant_Effect == combs$Type2[i],Pop==pop) %>% pull(RXY)
    
    # Calculate differences
    diff=Type2-Type1
    # Now calculate statistics for testing
    diff_mean=mean(diff)
    diff_se=jackknife_se(diff)
    
    # Compute paired t-statistic
    t_stat <- diff_mean / diff_se
    
    # One-tailed p-value (expecting mean_diff > 0)
    p_value <- pt(t_stat, df = n_blocks - 1)  
    
    print(paste0(pop,": ",combs$Type2[i]," < ",combs$Type1[i],
                 ": t = ",round(t_stat,3),
                 ", P = ",round(p_value,3)))
    
  } 
}


############################################
###          Mutations in ROH            ###
############################################

# Again, we do the same thing for RXY to see if difference mutation classes are
# more or less likely to fall inside of ROH.

# Refer to code above for comments on block definition code.

# LHISI

tmp <- effects %>% 
  filter(fate_lhisi == "segregating") %>% 
  mutate(ROH_50 = ifelse(Coverage_lhisi > 0.5,"Upper","Lower"))

# Define blocks
tmp <- tmp[order(tmp$Scaffold, tmp$Position), ]
tmp$Block <- NA
scafs <- unique(tmp$Scaffold)
total_blocks <- 93
total_entries <- nrow(tmp)
# Calculate the rough number of entries per block
rough_block_size <- ceiling(total_entries / total_blocks)
for(i in 1:length(scafs)){
  entries <- tmp[tmp$Scaffold==scafs[i],]
  blocks_in <- ceiling(nrow(entries)/rough_block_size)
  blocks_assign <- sort(rep(1:blocks_in, length.out=nrow(entries)))
  if(i == 1){
    
    tmp$Block[tmp$Scaffold==scafs[i]] <- blocks_assign
    
  } else {
    
    blocks_assign <- blocks_assign + max(tmp$Block[tmp$Scaffold==scafs[i-1]])
    tmp$Block[tmp$Scaffold==scafs[i]] <- blocks_assign
    
  }
}

n_blocks <- max(tmp$Block)

# Jackknife, blocks
for(i in 1:n_blocks){
  
  # Remove block, calculate proportion of loci per type in low ROH regions
  assign(paste0("jack_roh_",i),
         tmp %>% filter(Block != i) %>%
           group_by(Variant_Effect) %>%
           dplyr::summarise(Ratio=sum(ROH_50=="Lower")/n(),Jack=i))

}

# Collate results 
table_names <- ls(pattern = "^jack_roh_")
table_list <- mget(table_names)
roh_lhisi <- do.call(rbind, table_list) %>% mutate(Pop="LHISI")
rm(list=ls(pattern = "^jack_roh_"))

# LHIP

tmp <- effects %>% 
  filter(fate_lhip == "segregating") %>% 
  mutate(ROH_50 = ifelse(Coverage_lhip > 0.5,"Upper","Lower"))

# Define blocks
tmp <- tmp[order(tmp$Scaffold, tmp$Position), ]
tmp$Block <- NA
scafs <- unique(tmp$Scaffold)
total_blocks <- 92
total_entries <- nrow(tmp)
rough_block_size <- ceiling(total_entries / total_blocks)
for(i in 1:length(scafs)){
  entries <- tmp[tmp$Scaffold==scafs[i],]
  blocks_in <- ceiling(nrow(entries)/rough_block_size)
  blocks_assign <- sort(rep(1:blocks_in, length.out=nrow(entries)))
  if(i == 1){
    
    tmp$Block[tmp$Scaffold==scafs[i]] <- blocks_assign
    
  } else {
    
    blocks_assign <- blocks_assign + max(tmp$Block[tmp$Scaffold==scafs[i-1]])
    tmp$Block[tmp$Scaffold==scafs[i]] <- blocks_assign
    
  }
}

n_blocks <- max(tmp$Block)

# Jackknife, blocks
for(i in 1:n_blocks){
  
  # Remove block, calculate proportion of loci per type in low ROH regions
  assign(paste0("jack_roh_",i),
         tmp %>% filter(Block != i) %>%
           group_by(Variant_Effect) %>%
           dplyr::summarise(Ratio=sum(ROH_50=="Lower")/n(),Jack=i))
  
}

# Collate results 
table_names <- ls(pattern = "^jack_roh_")
table_list <- mget(table_names)
roh_lhip <- do.call(rbind, table_list) %>% mutate(Pop="LHIP")
rm(list=ls(pattern = "^jack_roh_"))

# And combine pops
jackknifes <- rbind(roh_lhisi,
                    roh_lhip)
jackknifes$Variant_Effect <- factor(jackknifes$Variant_Effect,
                                    levels=c("HIGH","MODERATE","LOW","MODIFIER"),
                                    labels=c("HIGH","MODERATE","LOW","MODIFIER"))

# Plot the distributions of Rxy for each
jackknife_se <- function(x){
  sqrt((length(x) - 1) / length(x) * sum((x - mean(x))^2))
}

plotting_results <- jackknifes %>%
  group_by(Variant_Effect,Pop) %>%
  dplyr::summarise(Mean=mean(Ratio),
                   SE=jackknife_se(Ratio),
                   Upper95=Mean+1.96*SE,
                   Lower95=Mean-1.96*SE) %>%
  ungroup

p <- plotting_results %>% 
  ggplot() +
  geom_errorbar(aes(y=Pop,xmax=Mean+SE,xmin=Mean-SE,group=Variant_Effect),
                size=0.75,position=position_dodge(width=0.6),width=0,
                filter(plotting_results,Pop=="LHISI")) + 
  geom_point(aes(y=Pop,x=Mean,fill=Variant_Effect,group=Variant_Effect),
             position=position_dodge(width=0.6),pch=21,size=4.5,
             filter(plotting_results,Pop=="LHISI"))  +
  scale_fill_manual(values=rev(new_colours[,3])) + 
  new_scale_fill() +
  geom_errorbar(aes(y=Pop,xmax=Mean+SE,xmin=Mean-SE,group=Variant_Effect),
                size=0.75,position=position_dodge(width=0.6),width=0,
                filter(plotting_results,Pop=="LHIP")) + 
  geom_point(aes(y=Pop,x=Mean,fill=Variant_Effect),
             position=position_dodge(width=0.6),pch=21,size=4.5,
             filter(plotting_results,Pop=="LHIP")) + 
  scale_fill_manual(values=rev(new_colours[,2])) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        legend.position="none"
  ) + 
  scale_x_continuous(
    limits=c(0.49,0.78),
    breaks=c(0.55,0.65,0.75),
    labels=c(55,65,75)) +
  labs(x="% variants in\nlow ROH regions") ; p

# Plot

png(paste0(FIGURE_DIR,"/muts_in_roh_",format(Sys.time(),"%Y%m%d"),".png"),
    width=6,height=3,units='in',res=300)
plot(p)
dev.off()

# We do the same paired t-test approach we did above, except this time we
# expect the more impact ful mutation classes to have a greater value than the
# one below, so we calculate the difference in the other direction.

# Generate a data.frame of combinations manually then loop over combinations.

combs <- data.frame(Type1=c("MODIFIER","MODIFIER","MODIFIER",
                            "LOW","LOW",
                            "MODERATE"),
                    Type2=c("LOW","MODERATE","HIGH",
                            "MODERATE","HIGH",
                            "HIGH")) %>% arrange(Type2)
pops <- c("LHISI","LHIP")

for(pop in pops){
  for(i in 1:nrow(combs)){
    
    # Get vectors of values
    Type1=jackknifes %>% arrange(Jack) %>% filter(Variant_Effect == combs$Type1[i],Pop==pop) %>% pull(Ratio)
    Type2=jackknifes %>% arrange(Jack) %>% filter(Variant_Effect == combs$Type2[i],Pop==pop) %>% pull(Ratio)
    
    # Calculate differences
    diff=Type2-Type1
    # Now calculate statistics for testing
    diff_mean=mean(diff)
    diff_se=jackknife_se(diff)
    
    # Compute paired t-statistic
    t_stat <- diff_mean / diff_se
    
    # One-tailed p-value (expecting mean_diff > 0)
    p_value <- pt(t_stat, df = n_blocks - 1,lower.tail=F)  
    
    print(paste0(pop,": ",combs$Type2[i]," > ",combs$Type1[i],
                 ": t = ",round(t_stat,3),
                 ", P = ",round(p_value,3)))
    
  } 
}
