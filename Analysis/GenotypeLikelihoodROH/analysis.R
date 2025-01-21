### Script to analyse and plot ROH data for lcWGS

# Libraries

library(dplyr)
library(magrittr)
library(ggplot2)
library(RZooRoH)
library(tidyr)
library(patchwork)
library(cowplot)

# Set environment

source("../../config.R")
source(paste0(DATA_DIR,"/colours.R"))
setwd(paste0(WORKING_DIR,"/GenotypeLikelihoodROH"))

# My favourite function

`%ni%` <- Negate(`%in%`)

# Results

load("all_results/all_model_results.RData")

# Metadata

pop="all"
names <- readLines(paste0("../../",pop,"_bam_list.txt"))
names <- data.frame(Sample=gsub("\\.bam","",gsub(".*/","",names)),
                    id=1:length(names))
details <- merge(names,read.table(paste0(REF_DIR,"/wgs_sample_details.txt"),
                                  header=T,stringsAsFactors=F,sep="\t"))
details$Pop <- factor(details$Pop,levels=c("Wild","LHIP","LHISI"))
details <- details[order(details$Pop),]
details$Plot_Order <- 1:nrow(details)
lengths <- read.table(paste0(REF_DIR,"/autosome_regions.bed"))[,c(1,3)]
colnames(lengths) <- c("Scaffold","total_length")
lengths$chrom <- 1:16

# Select model of K=9,rate=2

results <- all_K9_rate2_results
d <- results@hbdseg

# Add in metadata

d <- merge(d,details,by='id')
d <- merge(d,lengths,by='chrom')

if(!file.exists("figures")){
  dir.create("figures")
}

# Make figures of ROH across genome

for(i in 1:16){
  
  d_i <- d[d$chrom == i,]
  max_length <- max(d_i$total_length)
  
  p <- ggplot(d_i) + 
    geom_rect(aes(xmin=start_pos,xmax=end_pos,ymin=Plot_Order-0.3,ymax=Plot_Order+0.3,fill=Pop)) + 
    scale_x_continuous(limits=c(0-1e6,max_length+1e6),expand=c(0,0),
                       breaks=seq(from=0,to=round(max(d_i$end_pos),-8),length.out=5)) + 
    geom_hline(aes(yintercept=Plot_Order,color=Pop)) +
    scale_y_continuous(breaks=c(1:max(d_i$id)),
                       labels=unique(d_i$Sample)[order(unique(d_i$Plot_Order))],
                       limits=c(0.4,max(d_i$Plot_Order)+0.6)) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks=element_blank(),
          #        axis.title.x=element_blank(),
          #        axis.text.x=element_blank()
    ) + 
    xlab("Position") + 
    ggtitle(d_i$Scaffold[1]) +
    my_fill_3 + my_col_3
  
  png(paste0("figures/all_roh_",d_i$Scaffold[1],".png"),res=300,width=6.5,height=5,units='in')
  plot(p)
  dev.off()
  
}

# Plot Fis

fis <- cbind(details[order(details$id),],1-results@realized[,9])
colnames(fis)[ncol(fis)] <- "F_ROH"

p <- ggplot(fis,aes(y=F_ROH,x=Pop)) + 
  stat_summary(aes(x=Pop,y=F_ROH),
               fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5) +
  geom_point(pch=21,
             size=4,
             position=position_jitterdodge(dodge.width=1,jitter.width=0.5),
             aes(fill=Pop)) +  
  scale_y_continuous(limits=c(0,1)) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks=element_blank(),
        panel.grid = element_blank()) + 
  my_fill_3 + 
  labs(y=expression(F[ROH]))

# Save figure

png(paste0(FIGURE_DIR,"/froh_",format(Sys.time(),"%Y%m%d"),".png"),
    width=4.5,height=4.5,res=300,units='in')
plot(p)
dev.off()

# Test FROH difference

# We test excluding/including C01220 outlier

fis %>% filter(Pop != "Wild") %>% t.test(F_ROH~Pop,data=.)
fis %>% filter(Pop != "Wild",Sample != "C01220") %>% t.test(F_ROH~Pop,data=.)

# Plotting CM056998.1 as an example of ROH locations

d_i <- d[d$Scaffold == "CM056998.1",]
max_length <- lengths$V3[lengths$V1 == "CM056998.1"]

p <- ggplot(d_i) + 
  geom_rect(aes(xmin=start_pos,xmax=end_pos,ymin=Plot_Order-0.3,ymax=Plot_Order+0.3,fill=Pop)) + 
  scale_x_continuous(
#    limits=c(0-1e6,max_length+1e6),expand=c(0,0),
    breaks=seq(from=0,to=round(max(d_i$end_pos),-8),length.out=5),
    labels=c(expression(0),
                              expression(0.5%*%10^8),
                              expression(1.0%*%10^8),
                              expression(1.5%*%10^8),
                              expression(2.0%*%10^8))) + 
  geom_hline(aes(yintercept=Plot_Order,color=Pop)) +
  scale_y_continuous(breaks=c(1:max(d_i$id)),
                     labels=unique(d_i$Sample)[order(unique(d_i$Plot_Order))],
                     limits=c(0.4,max(d_i$Plot_Order)+0.6)) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        #        axis.title.x=element_blank(),
        #        axis.text.x=element_blank(),
        legend.position="none"
  ) + 
  xlab("Position") + 
  my_fill_3 + my_col_3 + 
  annotate(geom="rect",
           xmin=1e8,xmax=1.18e8,
           ymin=min(d_i$Plot_Order)-0.5,ymax=max(d_i$Plot_Order)+0.5,
           colour="black", linetype="dashed" ,fill=NA)

# Save figure

png(paste0(FIGURE_DIR,"/roh_location_example_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=7,height=2.5,units='in')
    plot(p)
dev.off()

# Plot length and number of ROH

# Create length classes

d$length_class <- cut(d$length,
                      breaks=c(0, 100000, 300000, 500000, 1000000, 1000000000),
                      labels=c('0 - 0.1', '0.1 - 0.3', '0.3 - 0.5', '0.5 - 1','>1'))

d$length_class <- factor(d$length_class,
                         levels=c('0 - 0.1', '0.1 - 0.3', '0.3 - 0.5', '0.5 - 1','>1'),
                         labels=c('0 - 0.1', '0.1 - 0.3', '0.3 - 0.5', '0.5 - 1','>1'))


p1 <- d %>%
  group_by(length_class,Sample) %>%
  dplyr::summarise(total_roh = sum(length),
                   Pop=first(Pop),
                   Metric="Length") %>%
  ggplot() +
  stat_summary(aes(x=length_class,y=total_roh,group=Pop),
               fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 1,position=position_dodge()) +
  geom_point(pch=21,size=3.5,
             position=position_jitterdodge(jitter.width=0.2,dodge.width=1),
             aes(x=length_class,y=total_roh,fill=Pop,group=Pop)) +
  scale_y_continuous(
    trans="log10",
    limits=c(7e5,1.5e9),
    breaks=c(1e6,1e7,1e8,1e9),
    labels=c(expression(1%*%10^6),
             expression(1%*%10^7),
             expression(1%*%10^8),
             expression(1%*%10^9))
  ) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        panel.border=element_rect(colour="black",fill=NA),
        axis.line = element_line(),
        axis.text.x=element_text(angle=30,hjust=1,vjust=1),
        panel.grid=element_blank()) +
  ylab("Length of genome\nin ROH (bp)") + xlab("Length (Mbp)") +
  my_fill_3

p2 <- d %>%
  group_by(length_class,Sample) %>%
  dplyr::summarise(n_roh = n(),
                   Pop=first(Pop),
                   Metric="Number") %>%
  ggplot() +
  stat_summary(aes(x=length_class,y=n_roh,group=Pop),
               fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 1,position=position_dodge()) +
  geom_point(pch=21,size=3.5,
             position=position_jitterdodge(jitter.width=0.2,dodge.width = 1),
             aes(x=length_class,y=n_roh,fill=Pop,group=Pop)) +
  scale_y_continuous(
  ) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        panel.border=element_rect(colour="black",fill=NA),
        axis.line = element_line(),
        axis.text.x=element_text(angle=30,hjust=1,vjust=1),
        panel.grid=element_blank()) +
  ylab("Number of ROH") + xlab("Length (Mbp)") +
  my_fill_3

p <- (p1 + theme(axis.title.x=element_blank(),axis.text.x=element_blank())) /
  p2 +
  plot_layout(guides="collect")

# Save figure

png(paste0(FIGURE_DIR,"/roh_tally_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=7,height=6,units='in')
    plot(p)
dev.off()

# Now the tests

# Number without filtering C01220

d %>% filter(Pop != "Wild") %>%
  group_by(length_class,Sample) %>%
  dplyr::summarise(n_roh = n(),
                   Pop=first(Pop)) %>%
  dplyr::select(Pop,length_class,n_roh) %>%
  pivot_wider(names_from=Pop,values_from=n_roh) %>%
  group_by(length_class) %>%
  dplyr::mutate(lhisi_mean = mean(unlist(LHISI)),
                lhip_mean = mean(unlist(LHIP)),
                p_value = t.test(unlist(LHISI), unlist(LHIP))$p.value,
                t_value = t.test(unlist(LHISI), unlist(LHIP))$statistic,
                df = t.test(unlist(LHISI), unlist(LHIP))$parameter
  )

# Filter C01220

d %>% filter(Pop != "Wild", Sample != "C01220") %>%
  group_by(length_class,Sample) %>%
  dplyr::summarise(n_roh = n(),
                   Pop=first(Pop)) %>%
  dplyr::select(Pop,length_class,n_roh) %>%
  pivot_wider(names_from=Pop,values_from=n_roh) %>%
  group_by(length_class) %>%
  dplyr::mutate(lhisi_mean = mean(unlist(LHISI)),
                lhip_mean = mean(unlist(LHIP)),
                p_value = t.test(unlist(LHISI), unlist(LHIP))$p.value,
                t_value = t.test(unlist(LHISI), unlist(LHIP))$statistic,
                df = t.test(unlist(LHISI), unlist(LHIP))$parameter
  )

# Log length_sum without filtering C01220

d %>% filter(Pop != "Wild") %>%
  group_by(length_class,Sample) %>%
  dplyr::summarise(length_sum = log(sum(length)),
                   Pop=first(Pop)) %>%
  dplyr::select(Pop,length_class,length_sum) %>%
  pivot_wider(names_from=Pop,values_from=length_sum) %>%
  group_by(length_class) %>%
  dplyr::mutate(lhisi_mean = mean(unlist(LHISI)),
                lhip_mean = mean(unlist(LHIP)),
                p_value = t.test(unlist(LHISI), unlist(LHIP))$p.value,
                t_value = t.test(unlist(LHISI), unlist(LHIP))$statistic,
                df = t.test(unlist(LHISI), unlist(LHIP))$parameter
  )

# Filter C01220

d %>% filter(Pop != "Wild", Sample != "C01220") %>%
  group_by(length_class,Sample) %>%
  dplyr::summarise(length_sum = log(sum(length)),
                   Pop=first(Pop)) %>%
  dplyr::select(Pop,length_class,length_sum) %>%
  pivot_wider(names_from=Pop,values_from=length_sum) %>%
  group_by(length_class) %>%
  dplyr::mutate(lhisi_mean = mean(unlist(LHISI)),
                lhip_mean = mean(unlist(LHIP)),
                p_value = t.test(unlist(LHISI), unlist(LHIP))$p.value,
                t_value = t.test(unlist(LHISI), unlist(LHIP))$statistic,
                df = t.test(unlist(LHISI), unlist(LHIP))$parameter
  )

# Finally, save per individual ROH bed files for use later

samples <- unique(d$Sample)

for(i in 1:length(samples)){
  
  tmp <- d[d$Sample == samples[i],]
  write.table(data.frame(scaffold=tmp$Scaffold,
                         start=tmp$start_pos-1,
                         end=tmp$end_pos-1),
              file=paste0(samples[i],"_rohs.bed"),col.names=F,row.names=F,sep="\t",quote=F)
  
  
}

# Use LHISI Fis values to estimate Ne in captivity
# Slightly hacky solution

# Function to compute change in F over 16 generations as a function of Ne
f1 <- function(Ne){
  1-((1-(1/(2*Ne)))^16)
}

# Range to evaluate function over
range=seq(1,100,0.0001)

# Fis from LHISI
fis_lhisi <- fis %>% filter(Pop=="LHISI") %>% pull(F_ROH)

# Wild Fis estimate, mean of two wild individuals
fis_wild <- fis %>% filter(Pop=="Wild") %>% pull(F_ROH) %>% mean

# Storage vector
vec <- c()

# Bootstrap to get confidence intervals
for(i in 1:10000){
  
  if(i %% 100 == 0){print(i)}
  
  #Resample
  resampled_mean <- sample(fis_wild,size=length(fis_wild),replace=T)
  
  # What is the difference between jackknifed LHISI F and measured wild F
  diff <- abs(fis_wild - mean(fis_lhisi[-i]))
  
  # Which Ne results in the same F difference after 16 generations?
  vec[i] <- range[which.min(abs(diff-f1(Ne=range)))]
  
}

result <- data.frame(Ne=vec) %>%
  summarise(Mean=mean(Ne),
            SE=sqrt((n() - 1) / n() * sum((Ne - mean(Ne))^2)),
            Upper95=Mean+1.96*SE,
            Lower95=Mean-1.96*SE) ; result

# Now examine the bias generated by not using population specific allele frequencies

# Read in all data objects

all <- load("all_results/all_model_results.RData")
lhip <- load("lhip_results/lhip_model_results.RData")
lhisi <- load("lhisi_results/lhisi_model_results.RData")

# Compare parameters per population on F

# Make vector of model parameters explored

rates=c(2,5,10)
classes=c(6:10)

# Empty data.frame

d <- data.frame()

# Loop over results objects, extracting FROH and adding columns for model parameters

for(rate in 1:length(rates)){
  
  for(class in 1:length(classes)){
    
    tmp <- eval(parse(text=paste0("all_K",classes[class],"_rate",rates[rate],"_results")))
    tmp <- cbind(details[order(details$id),],1-tmp@realized[,classes[class]])
    colnames(tmp)[ncol(tmp)] <- "F_ROH"
    tmp$K = classes[class] ; tmp$Rate <- rates[rate]
    d <- rbind(d,tmp)
    
  }
  
}

# Plot

labels <- c(
  `2` = "r = 2",
  `5` = "r = 5",
  `10` = "r = 10"
)

p <- ggplot(d,aes(x=K,y=F_ROH,group=Sample,shape=Sex,colour=Pop)) + 
  geom_path() + 
  geom_point() + 
  facet_wrap(~Rate,labeller=as_labeller(labels)) + 
  theme_bw() + 
  theme(axis.ticks=element_blank()) + 
  labs(y=expression(F[ROH])) + 
  my_col_3

# Save figure

png(paste0(FIGURE_DIR,"/roh_parameter_comp_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=7,height=6,units='in')
plot(p)
dev.off()

# FROH from different base populations

# Empty data.frame

output <- data.frame(fs=numeric(),names=character(),
                     k=numeric(),rate=numeric(),pop=character())

pops <- c("all","lhisi","lhip")

for(i in 1:length(pops)){
  
  fs <- eval(parse(text=paste0("1- ",pops[i],"_K9_rate2_results@realized[,9]")))
  names <- readLines(paste0(paste0("../../",pops[i],"_bam_list.txt")))
  names <- data.frame(Sample=gsub("\\.bam","",gsub(".*/","",names)),
                      id=1:length(names))
  
  output <- rbind(output,
                  data.frame(fs=fs,names=names$Sample,pop=pops[i]))
  
}

# We know from above that all models are very similar
# For all comparisons, we are just looking at K=9, rate=2

# Calculate a mean difference in f per individual and plot

p <- output %>% filter(names != "VAN2" & names != "PAUL4") %>%
  pivot_wider(names_from=pop,values_from=fs) %>%
  mutate(LHISI=all-lhisi,LHIP=all-lhip) %>%
  dplyr::select(-all,-lhisi,-lhip) %>%
  gather(-names,key="population",value="Difference") %>%
  drop_na() %>%
  ggplot(aes(x=population,y=Difference,fill=population)) + 
  stat_summary(aes(x=population,y=Difference),
               fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5) +
  geom_point(position=position_jitter(width=0.15),pch=21,size=4) + 
  theme_bw() + 
  my_fill_2 +
  labs(x="Population") +
  theme(axis.ticks=element_blank(),
        title=element_text(size=9),
        legend.position="none") + 
  scale_y_continuous(limits=c(-0.1,0.1)) + 
  scale_x_discrete(labels=c("Hybrid","Captive"))

# Save figure

png(paste0(FIGURE_DIR,"/roh_parameter_comp_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=7,height=6,units='in')
plot(p)
dev.off()

# Now also plot FROH

d_lhisi <- data.frame(Pop="LHISI",
                      F_ROH=1-lhisi_K9_rate2_results@realized[,9])
d_lhip <- data.frame(Pop="LHIP",
                     F_ROH=1-lhip_K9_rate2_results@realized[,9])
fis_specific <- rbind(d_lhisi,d_lhip)

p <- ggplot(fis_specific,aes(y=F_ROH,x=Pop)) + 
  stat_summary(aes(x=Pop,y=F_ROH),
               fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.5) +
  geom_point(pch=21,
             size=4,
             position=position_jitterdodge(dodge.width=1,jitter.width=0.5),
             aes(fill=Pop)) +  
  scale_y_continuous(limits=c(0,1)) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks=element_blank(),
        panel.grid = element_blank()) + 
  my_fill_3 + 
  labs(y=expression(F[ROH]))

png(paste0(FIGURE_DIR,"/froh_with_pop_specific",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=7,height=6,units='in')
plot(p)
dev.off()


# What is the difference in F between lhip
# and lhisi when only considering population
# specific allele frequencies

d <- output %>% filter(names != "VAN2" & names != "PAUL4" & pop != "all")

details <- read.table(paste0(REF_DIR,"/wgs_sample_details.txt"),
                      header=T,stringsAsFactors = F)
colnames(details) <- c("names","sex","population")
output <- merge(d,details)

# Does statistical significance disappear when using these population-specific models

output %>% t.test(fs~pop,data=.)
output %>% filter(names != "C01220") %>% t.test(fs~pop,data=.)

# No

# Now redo the tests for comparing length and number of ROH between hybrid and wild

d_lhisi <- lhisi_K9_rate2_results@hbdseg
d_lhisi$Pop <- "LHISI"
d_lhisi$id <- paste0("LHISI_",d_lhisi$id)

d_lhip <- lhip_K9_rate2_results@hbdseg
d_lhip$Pop <- "LHIP"
d_lhip$id <- paste0("LHIP_",d_lhip$id)

# FIND C01220
c01220 <- d_lhip %>% group_by(id) %>% summarise(froh=sum(length)) %>% filter(froh==max(froh)) %>% pull(id)

d <- rbind(d_lhisi,d_lhip)

d$length_class <- cut(d$length,
                      breaks=c(0, 100000, 300000, 500000, 1000000, 1000000000),
                      labels=c('0 - 0.1', '0.1 - 0.3', '0.3 - 0.5', '0.5 - 1','>1'))

d$length_class <- factor(d$length_class,
                         levels=c('0 - 0.1', '0.1 - 0.3', '0.3 - 0.5', '0.5 - 1','>1'),
                         labels=c('0 - 0.1', '0.1 - 0.3', '0.3 - 0.5', '0.5 - 1','>1'))

# Number without filtering C01220

d %>% 
  group_by(length_class,id) %>%
  dplyr::summarise(n_roh = n(),
                   Pop=first(Pop)) %>%
  dplyr::select(Pop,length_class,n_roh) %>%
  pivot_wider(names_from=Pop,values_from=n_roh) %>%
  group_by(length_class) %>%
  dplyr::mutate(lhisi_mean = mean(unlist(LHISI)),
                lhip_mean = mean(unlist(LHIP)),
                p_value = t.test(unlist(LHISI), unlist(LHIP))$p.value,
                t_value = t.test(unlist(LHISI), unlist(LHIP))$statistic,
                df = t.test(unlist(LHISI), unlist(LHIP))$parameter
  )

# Filter C01220

d %>% filter(id != c01220) %>%
  group_by(length_class,id) %>%
  dplyr::summarise(n_roh = n(),
                   Pop=first(Pop)) %>%
  dplyr::select(Pop,length_class,n_roh) %>%
  pivot_wider(names_from=Pop,values_from=n_roh) %>%
  group_by(length_class) %>%
  dplyr::mutate(lhisi_mean = mean(unlist(LHISI)),
                lhip_mean = mean(unlist(LHIP)),
                p_value = t.test(unlist(LHISI), unlist(LHIP))$p.value,
                t_value = t.test(unlist(LHISI), unlist(LHIP))$statistic,
                df = t.test(unlist(LHISI), unlist(LHIP))$parameter
  )

# Log length_sum without filtering C01220

d %>% 
  group_by(length_class,id) %>%
  dplyr::summarise(length_sum = log(sum(length)),
                   Pop=first(Pop)) %>%
  dplyr::select(Pop,length_class,length_sum) %>%
  pivot_wider(names_from=Pop,values_from=length_sum) %>%
  group_by(length_class) %>%
  dplyr::mutate(lhisi_mean = mean(unlist(LHISI)),
                lhip_mean = mean(unlist(LHIP)),
                p_value = t.test(unlist(LHISI), unlist(LHIP))$p.value,
                t_value = t.test(unlist(LHISI), unlist(LHIP))$statistic,
                df = t.test(unlist(LHISI), unlist(LHIP))$parameter
  )

# Filter C01220

d %>% filter(id != c01220) %>%
  group_by(length_class,id) %>%
  dplyr::summarise(length_sum = log(sum(length)),
                   Pop=first(Pop)) %>%
  dplyr::select(Pop,length_class,length_sum) %>%
  pivot_wider(names_from=Pop,values_from=length_sum) %>%
  group_by(length_class) %>%
  dplyr::mutate(lhisi_mean = mean(unlist(LHISI)),
                lhip_mean = mean(unlist(LHIP)),
                p_value = t.test(unlist(LHISI), unlist(LHIP))$p.value,
                t_value = t.test(unlist(LHISI), unlist(LHIP))$statistic,
                df = t.test(unlist(LHISI), unlist(LHIP))$parameter
  )

# Also redo this plot

p1 <- d %>%
  group_by(length_class,id) %>%
  dplyr::summarise(total_roh = sum(length),
                   Pop=first(Pop),
                   Metric="Length") %>%
  ggplot() +
  stat_summary(aes(x=length_class,y=total_roh,group=Pop),
               fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 1,position=position_dodge()) +
  geom_point(pch=21,size=3.5,
             position=position_jitterdodge(jitter.width=0.2,dodge.width=1),
             aes(x=length_class,y=total_roh,fill=Pop,group=Pop)) +
  scale_y_continuous(
    trans="log10",
    limits=c(7e5,1.5e9),
    breaks=c(1e6,1e7,1e8,1e9),
    labels=c(expression(1%*%10^6),
             expression(1%*%10^7),
             expression(1%*%10^8),
             expression(1%*%10^9))
  ) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        panel.border=element_rect(colour="black",fill=NA),
        axis.line = element_line(),
        axis.text.x=element_text(angle=30,hjust=1,vjust=1),
        panel.grid=element_blank()) +
  ylab("Length of genome\nin ROH (bp)") + xlab("Length (Mbp)") +
  my_fill_3


p2 <- d %>%
  group_by(length_class,id) %>%
  dplyr::summarise(n_roh = n(),
                   Pop=first(Pop),
                   Metric="Number") %>%
  ggplot() +
  stat_summary(aes(x=length_class,y=n_roh,group=Pop),
               fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 1,position=position_dodge()) +
  geom_point(pch=21,size=3.5,
             position=position_jitterdodge(jitter.width=0.2,dodge.width = 1),
             aes(x=length_class,y=n_roh,fill=Pop,group=Pop)) +
  scale_y_continuous(
  ) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        panel.border=element_rect(colour="black",fill=NA),
        axis.line = element_line(),
        axis.text.x=element_text(angle=30,hjust=1,vjust=1),
        panel.grid=element_blank()) +
  ylab("Number of ROH") + xlab("Length (Mbp)") +
  my_fill_3

p <- (p1 + theme(axis.title.x=element_blank(),axis.text.x=element_blank())) /
  p2 +
  plot_layout(guides="collect")

# Save figure

png(paste0(FIGURE_DIR,"/roh_tally_different_base_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=7,height=6,units='in')
plot(p)
dev.off()

# Calculate whether mean depth and Froh are correlated

depth <- read.table("../AnalysingWGSCoverage/coverage_counts.txt", header=F,stringsAsFactors = F,sep="\t")
depth <- depth[grep("CM",depth$V1),]
colnames(depth) <- c("Chrom","Start","End","Mean","Sample")
plot_data <- depth %>%
  group_by(Sample) %>%
  dplyr::summarise(Mean_Coverage = mean(Mean)) %>%
  merge(.,fis,by="Sample")
p <- plot_data %>%
  ggplot(aes(x=Mean_Coverage,y=F_ROH)) + 
  geom_smooth(method="lm",colour="black") +
  geom_point(size=3,aes(shape=Sex,fill=Pop)) + 
  scale_shape_manual(values=c(21,24)) +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill="lightgrey"))) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        panel.grid = element_blank()) + 
  my_fill_3 + 
  
  labs(y=expression(F[ROH]),x="Mean sample depth\nof coverage")

png(paste0(FIGURE_DIR,"/depth_and_froh_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=5,height=5,units='in')
plot(p)
dev.off()

# Also test

summary(lm(F_ROH~Mean_Coverage,data=plot_data))
