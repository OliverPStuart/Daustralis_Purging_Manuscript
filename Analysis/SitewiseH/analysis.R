
### Script to analyse and plot type-specific H from lcWGS

setwd("/Volumes/Alter/Daus_WGS_Paper/Analysis/SitewiseH/")

# Libraries

library(dplyr)
library(magrittr)
library(ggplot2)
library(ggnewscale)

# Set environment

source("../../config.R")
source(paste0(DATA_DIR,"/colours.R"))
setwd(paste0(WORKING_DIR,"/SitewiseH"))

# My favourite function

`%ni%` <- Negate(`%in%`)

# Colours
new_colours <- read.delim(paste0(REF_DIR,"/two_pop_2d_colour_scale.txt"),
                          header=F,stringsAsFactors = F)
colnames(new_colours) <- c("Effect","Captive","Hybrid","Wild")


### Part 1: Per-sample type-specific H
# We analyse each individual heterozygosity at different site types

# Data

data <- read.table("site_hets.txt",header=F,stringsAsFactors=F)[,1:5]
colnames(data) <- c("Sample","Depth","Type","HOM","HET")
data$H <- data$HET / (data$HOM + data$HET)
data <- data %>% dplyr::select(Sample,Depth,Type,H)
pops <- read.table("../../References/wgs_sample_details.txt",header=T,stringsAsFactors=F)

# Calculate H relative to MODIFIER, our neutral stand-in

data %>% 
  filter(Sample %ni% c("C01220")
         ) %>%
  mutate(Index=ifelse(Type=="MODIFIER",1,2),
         Type = factor(Type,levels=c("MODIFIER","LOW","MODERATE","HIGH"))) %>%
  arrange(Index) %>%
  group_by(Depth,Sample) %>%
  mutate(H_Rel = H/first(H)) %>%
  ungroup %>%
  filter(Type != "MODIFIER") %>%
  merge(.,pops) %>%
  ggplot(aes(x=Type,y=H_Rel)) + 
  stat_summary(aes(x=Type,y=H_Rel, group=Type),
               fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5) +
  geom_point(aes(x=Type,y=H_Rel, fill=Pop),pch=21,colour="black",
             size=4,position=position_jitterdodge(dodge.width=0.5,jitter.width=0.5)) +
  facet_grid(Depth~Pop) + 
  geom_hline(yintercept=1,linetype="dashed") + 
  my_fill_3 + 
  theme_bw() + 
  theme(axis.ticks=element_blank())

data %>% 
  filter(Sample %ni% c("C01220")
         ) %>%
  merge(.,pops) %>%
  ggplot(aes(x=Type,y=H)) + 
  stat_summary(aes(x=Type,y=H, group=Type),
               fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5) +
  geom_point(aes(x=Type,y=H, fill=Pop),pch=21,colour="black",
             size=4,position=position_jitterdodge(dodge.width=0.5,jitter.width=0.5)) +
  facet_grid(Depth~Pop) + 
  my_fill_3 + 
  theme_bw() + 
  theme(axis.ticks=element_blank())

# Well, let's test of these are different from 0
# Reshape to get Hs in different columns
data_reshaped <- data %>% 
  filter(Sample %ni% c("C01220","VAN2","PAUL4")) %>%
  tidyr::pivot_wider(names_from="Type",values_from="H") %>% 
  merge(.,pops)

# Test each type in each pop in each depth in turn
data_reshaped %>%
  group_by(Pop,Depth) %>%
  dplyr::mutate(MODIFIER_mean = mean(unlist(MODIFIER)),
                TEST_mean = mean(unlist(LOW)),
                p_value = t.test(unlist(LOW), unlist(MODIFIER),alternative="greater",paired=T)$p.value,
                t_value = t.test(unlist(LOW), unlist(MODIFIER),alternative="greater",paired=T)$statistic,
                df = t.test(unlist(LOW), unlist(MODIFIER),alternative="greater",paired=T)$parameter
  ) %>% select(MODIFIER_mean,TEST_mean,p_value,t_value,df) %>% distinct()

data_reshaped %>%
  group_by(Pop,Depth) %>%
  dplyr::mutate(MODIFIER_mean = mean(unlist(MODIFIER)),
                TEST_mean = mean(unlist(MODERATE)),
                p_value = t.test(unlist(MODERATE), unlist(MODIFIER),alternative="greater",paired=T)$p.value,
                t_value = t.test(unlist(MODERATE), unlist(MODIFIER),alternative="greater",paired=T)$statistic,
                df = t.test(unlist(MODERATE), unlist(MODIFIER),alternative="greater",paired=T)$parameter
  ) %>% select(MODIFIER_mean,TEST_mean,p_value,t_value,df) %>% distinct()

data_reshaped %>%
  group_by(Pop,Depth) %>%
  dplyr::mutate(MODIFIER_mean = mean(unlist(MODIFIER)),
                TEST_mean = mean(unlist(HIGH)),
                p_value = t.test(unlist(HIGH), unlist(MODIFIER),alternative="greater",paired=T)$p.value,
                t_value = t.test(unlist(HIGH), unlist(MODIFIER),alternative="greater",paired=T)$statistic,
                df = t.test(unlist(HIGH), unlist(MODIFIER),alternative="greater",paired=T)$parameter
  ) %>% select(MODIFIER_mean,TEST_mean,p_value,t_value,df) %>% distinct()

# Should we be comparing to the wild somehow???

data %>% 
  filter(Sample %ni% c("C01220")
  ) %>%
  merge(.,pops) %>%
  ggplot(aes(x=Pop,y=H)) + 
  stat_summary(aes(x=Pop,y=H, group=Type),
               fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5) +
  geom_point(aes(x=Pop,y=H, fill=Pop),pch=21,colour="black",
             size=4,position=position_jitterdodge(dodge.width=0.5,jitter.width=0.5)) +
  facet_grid(Depth~Type) + 
  my_fill_3 + 
  theme_bw() + 
  theme(axis.ticks=element_blank())

data %>% 
  filter(Sample %ni% c("C01220")
  ) %>%
  mutate(Index=ifelse(Type=="MODIFIER",1,2),
         Type = factor(Type,levels=c("MODIFIER","LOW","MODERATE","HIGH"))) %>%
  arrange(Index) %>%
  group_by(Depth,Sample) %>%
  mutate(H_Rel = H/first(H)) %>%
  ungroup %>%
  filter(Type != "MODIFIER") %>%
  merge(.,pops) %>%
  ggplot(aes(x=Pop,y=H_Rel)) + 
  stat_summary(aes(x=Pop,y=H_Rel, group=Type),
               fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5) +
  geom_point(aes(x=Pop,y=H_Rel, fill=Pop),pch=21,colour="black",
             size=4,position=position_jitterdodge(dodge.width=0.5,jitter.width=0.5)) +
  facet_grid(Depth~Type) + 
  geom_hline(yintercept=1,linetype="dashed") + 
  my_fill_3 + 
  theme_bw() + 
  theme(axis.ticks=element_blank())

### Part 2: Per-site type-specific H
# We analyse each subpopulations heterozygosity at different site types
# To make it congruent with previous results, we restrict it to segregating sites per population

hwe_lhisi <- read.table('lhisi_hwe.hwe.gz',header=T)[,c(1,2,10)] %>% filter(hetFreq > 0.01)
hwe_lhip <- read.table('lhip_hwe.hwe.gz',header=T)[,c(1,2,10)] %>% filter(hetFreq > 0.01)

high <- read.table("HIGH.bed")[,c(1,3)] ; high$Type <- "HIGH"
moderate <- read.table("MODERATE.bed")[,c(1,3)] ; moderate$Type <- "MODERATE"
low <- read.table("LOW.bed")[,c(1,3)] ; low$Type <- "LOW"
modifier <- read.table("MODIFIER.bed")[,c(1,3)] ; modifier$Type <- "MODIFIER"

sites <- rbind(high,moderate,low,modifier)
colnames(sites)[c(1,2)] <- c("Chromo","Position")

# How do we get mean heterozygosity for each site type?
# We use the same block jackknife approach as in the others

# LHISI

tmp <- merge(hwe_lhisi,sites,all.x=T)

# Order the data.frame by Scaffold and Position
tmp <- tmp[order(tmp$Chromo, tmp$Position), ]

# Initialise block column
tmp$Block <- NA

# Get scaffolds
scafs <- unique(tmp$Chromo)

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
  entries <- tmp[tmp$Chromo==scafs[i],]
  
  # How many blocks roughly in scaffold?
  blocks_in <- ceiling(nrow(entries)/rough_block_size)
  
  # Make a vector of block assignments
  blocks_assign <- sort(rep(1:blocks_in, length.out=nrow(entries)))
  
  # Add block assignments to data.frame
  # If it's not the first scaffold, add previous values to get correct block ID
  if(i == 1){
    
    tmp$Block[tmp$Chromo==scafs[i]] <- blocks_assign
    
  } else {
    
    blocks_assign <- blocks_assign + max(tmp$Block[tmp$Chromo==scafs[i-1]])
    tmp$Block[tmp$Chromo==scafs[i]] <- blocks_assign
    
  }
}

n_blocks <- max(tmp$Block)

# Jackknife, blocks
for(i in 1:n_blocks){
  
  # Remove block, calculate RXY, assign to temp data.frame
  assign(paste0("jack_ho_",i),
         tmp %>% filter(Block != i) %>% 
           group_by(Type) %>%
           dplyr::summarise(Ho=mean(hetFreq),Jack=i))
  
}

# Collate results 
table_names <- ls(pattern = "^jack_ho_")
table_list <- mget(table_names)
ho_lhisi <- do.call(rbind, table_list) %>% mutate(Pop="LHISI")
rm(list=ls(pattern = "^jack_ho_"))

# LHIP

tmp <- merge(hwe_lhip,sites,all.x=T)

# Order the data.frame by Scaffold and Position
tmp <- tmp[order(tmp$Chromo, tmp$Position), ]

# Initialise block column
tmp$Block <- NA

# Get scaffolds
scafs <- unique(tmp$Chromo)

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
  entries <- tmp[tmp$Chromo==scafs[i],]
  
  # How many blocks roughly in scaffold?
  blocks_in <- ceiling(nrow(entries)/rough_block_size)
  
  # Make a vector of block assignments
  blocks_assign <- sort(rep(1:blocks_in, length.out=nrow(entries)))
  
  # Add block assignments to data.frame
  # If it's not the first scaffold, add previous values to get correct block ID
  if(i == 1){
    
    tmp$Block[tmp$Chromo==scafs[i]] <- blocks_assign
    
  } else {
    
    blocks_assign <- blocks_assign + max(tmp$Block[tmp$Chromo==scafs[i-1]])
    tmp$Block[tmp$Chromo==scafs[i]] <- blocks_assign
    
  }
}

n_blocks <- max(tmp$Block)

# Jackknife, blocks
for(i in 1:n_blocks){
  
  # Remove block, calculate RXY, assign to temp data.frame
  assign(paste0("jack_ho_",i),
         tmp %>% filter(Block != i) %>% 
           group_by(Type) %>%
           dplyr::summarise(Ho=mean(hetFreq),Jack=i))
  
}

# Collate results 
table_names <- ls(pattern = "^jack_ho_")
table_list <- mget(table_names)
ho_lhip <- do.call(rbind, table_list) %>% mutate(Pop="LHIP")
rm(list=ls(pattern = "^jack_ho_"))

# And combine pops
jackknifes <- rbind(ho_lhisi,
                    ho_lhip)
jackknifes$Type <- factor(jackknifes$Type,
                                    levels=c("HIGH","MODERATE","LOW","MODIFIER"),
                                    labels=c("HIGH","MODERATE","LOW","MODIFIER"))

# Plot the distributions of Rxy for each
jackknife_se <- function(x){
  sqrt((length(x) - 1) / length(x) * sum((x - mean(x))^2))
}

plotting_results <- jackknifes %>%
  group_by(Type,Pop) %>%
  dplyr::summarise(Mean=mean(Ho),
                   SE=jackknife_se(Ho),
                   Upper95=Mean+1.96*SE,
                   Lower95=Mean-1.96*SE) %>%
  ungroup

p <- plotting_results %>%
  ggplot() + 
  geom_errorbar(aes(y=Pop,group=Type,xmax=Mean+SE,xmin=Mean-SE),
                size=1.25,position=position_dodge(width=0.6),width=0,
                filter(plotting_results,Pop=="LHISI")) + 
  geom_point(aes(y=Pop,group=Type,fill=Type,x=Mean),
             position=position_dodge(width=0.6), size=4.5,pch=21,
             filter(plotting_results,Pop=="LHISI")) + 
  scale_fill_manual(values=rev(new_colours[,3])) + 
  new_scale_fill() +
  geom_errorbar(aes(y=Pop,group=Type,xmax=Mean+SE,xmin=Mean-SE),
                size=1.25,position=position_dodge(width=0.6),width=0,
                filter(plotting_results,Pop=="LHIP")) + 
  geom_point(aes(y=Pop,group=Type,fill=Type,x=Mean),
             position=position_dodge(width=0.6), size=4.5,pch=21,
             filter(plotting_results,Pop=="LHIP")) + 
  scale_fill_manual(values=rev(new_colours[,2])) +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.text.x=element_text(size=10),
        #axis.text.y=element_blank(),
        #legend.position="none"
        ) + 
  scale_x_continuous(
    limits=c(0.41,0.51),
    breaks=seq(0.42,0.5,0.02)
    ) +
  labs(x=expression(H[O])) ; p  

# Now test!

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
    Type1=jackknifes %>% arrange(Jack) %>% filter(Type == combs$Type1[i],Pop==pop) %>% pull(Ho)
    Type2=jackknifes %>% arrange(Jack) %>% filter(Type == combs$Type2[i],Pop==pop) %>% pull(Ho)
    
    # Calculate differences
    diff=Type1-Type2
    # Now calculate statistics for testing
    diff_mean=mean(diff)
    diff_se=jackknife_se(diff)
    
    # Compute paired t-statistic
    t_stat <- diff_mean / diff_se
    
    # One-tailed p-value (expecting mean_diff > 0)
    p_value <- pt(t_stat, df = n_blocks - 1)  
    
    print(paste0(pop,": ",combs$Type2[i]," > ",combs$Type1[i],
                 ": t = ",round(t_stat,3),
                 ", P = ",round(p_value,3)))
    
  } 
}


