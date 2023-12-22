
### Script to analyse and plot depth in PAUL4 WGS sample

# Set environment

source("../../config.R")
setwd(paste0(WORKING_DIR,"/PAUL4"))

# Libraries

library(dplyr)
library(magrittr)
library(ggplot2)
library(Hmisc)
library(patchwork)
library(cowplot)

# Load data

H <- read.table("GENOME_WIDE_H.txt",header=F,sep="\t")
colnames(H) <- c("Chromosome","Start","End","Hom","Het")

# Calculate H and region coverage

H <- H %>%
  mutate(Coverage = (Hom+Het)/(End-Start),
         H=Het/(Hom+Het))

# First calculate whether the coverage of a region has anything to do with H

H %>%
  ggplot(aes(x=Coverage,y=H,fill=Chromosome)) + 
  geom_point(pch=21) + geom_smooth(method="lm") + 
  scale_fill_discrete(guide="none") + 
  facet_wrap(~Chromosome) + 
#  scale_x_continuous(limits=c(0.05,0.3),expand=c(0,0)) + 
  theme_bw()

# A little, so we're going to exclude regions below 0.1 coverage
# Now let's plot it all!

# Use scaffold length to get order of scaffolds, used for plotting
lengths <- H %>% select(Chromosome, End) %>%
  group_by(Chromosome) %>%
  filter(End==max(End)) %>%
  unique %>%
  ungroup() %>%
  mutate(Order=1:n()) %>%
  select(Chromosome,Order)

H <- merge(H,lengths)

H <- H[with(H, order(Chromosome, Start)),]

# Make bins for plotting
# Make continuous variable for plotting
H$Bin <- c(1:nrow(H))

# Get median and select corresponding bin number to get location of axis label
for(i in 1:length(unique(H$Chromosome))){
  
  if(i == 1){
    Locations <- c()
  }
  
  chrom=unique(H$Chromosome)[i]
  
  Locations <- c(Locations,
                 H %>% 
                   subset(Chromosome==chrom) %>%
                   summarise(Bin=median(Bin)) %>% unlist)
  
  if(i == length(unique(H$Chromosome))){
    Locations <- data.frame(Locations=Locations,Chromosome=unique(H$Chromosome))
  }
  
}

# Some formatting
H$Chromosome <- factor(H$Chromosome,levels=unique(H$Chromosome),labels=unique(H$Chromosome))

# Vector of colours
colours <- rep(c("springgreen4","lightgreen"),8)

# Set y-axis minimum and maximum
ymin=0
ymax=0.01

# Data.frame to make panel background rectangles
rect_breaks <- seq(from=ymin,to=ymax,by=(ymax-ymin)/5)
rect <- data.frame(ymin=rect_breaks[c(2,4)],
                   ymax=rect_breaks[c(3,5)],
                   xmax=Inf,xmin=-Inf)

# Make y_labels from this, max and min slightly nudged inward 
y_labels <- rect_breaks
y_labels[1] <- rect_breaks[1] + (rect_breaks[6]/100)
y_labels[6] <- rect_breaks[6] - (rect_breaks[6]/100)

# Make summary stats for reporting in the paper

weighted.se.mean <- function(x, w, na.rm = T){
  ## Remove NAs 
  if (na.rm) {
    i <- !is.na(x)
    w <- w[i]
    x <- x[i]
  }
  
  ## Calculate effective N and correction factor
  n_eff <- (sum(w))^2/(sum(w^2))
  correction = n_eff/(n_eff-1)
  
  ## Get weighted variance 
  numerator = sum(w*(x-weighted.mean(x,w))^2)
  denominator = sum(w)
  
  ## get weighted standard error of the mean 
  se_x = sqrt((correction * (numerator/denominator))/n_eff)
  return(se_x)
}

s_stats <- H %>%
  filter(Coverage > 0.1) %>%
  summarise(W_Mean = weighted.mean(H,Coverage),
            W_SE = weighted.se.mean(H,Coverage))

# Make a final label to get genome wide average

p1 <- H %>%
  filter(Coverage > 0.1) %>%
  ggplot(aes(x=Bin,group=Chromosome,y=H,colour=Chromosome,fill=Chromosome)) +
  theme_bw() +
  theme(axis.ticks=element_blank(),panel.grid=element_blank(),
        legend.position="none",
        axis.text.x=element_text(angle=35,hjust=1),
        axis.title.x=element_blank()) +
  geom_rect(inherit.aes=NULL,data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            colour="gray95",fill="gray95") +
  geom_line() +
  geom_area() +
  scale_x_continuous(expand=c(0,0),lim=c(1,nrow(H)),breaks=Locations$Locations,labels=Locations$Chromosome) +
  scale_y_continuous(breaks=y_labels,labels=rect_breaks) + 
  coord_cartesian(ylim=c(ymin,ymax),expand=F) + 
  scale_colour_manual(unique(H$Chromosome),values=colours) +
  scale_fill_manual(unique(H$Chromosome),values=colours) +
  labs(y="Heterozygosity",x="Chromosome") + 
  geom_hline(yintercept=s_stats$W_Mean,linetype="dashed") ; p1

png(paste0(FIGURE_DIR,"/genome_wide_H_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=10,height=5,units='in')
plot(p1)
dev.off()

# Now we're going to add to this the figure of genome_wide H of island species
# This is from published data on the Kakapo and the island foxes

islands <- read.table(paste0(DATA_DIR,"/20231222_genome_wide_H_islands.txt"),
                      header=T,stringsAsFactors = F,sep="\t")
islands <- rbind(islands,
      c("Dryococelus_australis","LHISI","Insecta",s_stats$W_Mean,"This study",NA))

p2 <- islands %>%
  ggplot(aes(x=as.numeric(Observed_H),fill=Common_Name)) + 
  geom_histogram(colour="black") + 
  scale_y_continuous(limits=c(0,7),expand=c(0,0)) +
  scale_fill_manual(values=c("grey94","grey70","springgreen3")) + 
  theme_bw() + 
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position=c(0.9,0.1),
        legend.title = element_blank(),
        legend.background = element_rect(colour="black"),
        panel.grid=element_blank()) + 
  labs(x="Genomic H"); p2

png(paste0(FIGURE_DIR,"/island_H_comparison_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=7,height=7,units='in')
plot(p2)
dev.off()
