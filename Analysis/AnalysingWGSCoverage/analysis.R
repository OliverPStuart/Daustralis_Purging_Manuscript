
### Script to analyse and plot depth in PAUL4 WGS sample

# Set environment

source("../../config.R")
setwd(paste0(WORKING_DIR,"/AnalysingWGSCoverage"))

# Libraries
library(dplyr)
library(magrittr)
library(ggplot2)

# Load data

depths <- read.table("PAUL4_WGS.regions.bed.gz")
colnames(depths) <- c("Scaffold","Start","End","Depth")

# Filter out minor scaffolds

depths <- depths[grep("CM",depths$Scaffold),]
unique(depths$Scaffold)

# Make categorical variable for sex chromosome

depths$Type <- ifelse(depths$Scaffold == "CM057009.1","X","Autosome")

# Make depth figure

p <- depths %>%
  ggplot(aes(x=Depth,fill=Type)) + 
  geom_histogram(colour="black") + 
  scale_x_continuous(limits=c(5,30),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,70000),expand=c(0,0)) +
  scale_fill_manual(values=c("#98d5fa","#2e6094")) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.grid=element_blank())

png(paste0(FIGURE_DIR,"/PAUL4_Depth.png"),res=300,width=5,height=3,units='in')
plot(p)
dev.off()

# Now analyse only autosomes
# How uniform is coverage? Can we use a universal filter?

grand_mean <- depths %>% 
  filter(Type == "Autosome") %>%
  summarise(Mean=mean(Depth)) %>%
  pull(Mean)

p <- depths %>% 
  filter(Type == "Autosome") %>%
  group_by(Scaffold) %>%
  summarise(Mean=mean(Depth),
            #SD=sd(Depth),
            Lower_95=quantile(Depth,probs=0.025),
            Upper_95=quantile(Depth,probs=0.975),
            Lower_99=quantile(Depth,probs=0.005),
            Upper_99=quantile(Depth,probs=0.995)) %>%
  ggplot() + 
  geom_vline(xintercept = grand_mean,linetype="dashed") +
  geom_errorbar(aes(y=Scaffold,x=Mean,xmin=Lower_95,xmax=Upper_95),
                width=0,size=1) + 
  geom_errorbar(aes(y=Scaffold,x=Mean,xmin=Lower_99,xmax=Upper_99),
                width=0,size=0.5) + 
  geom_point(aes(y=Scaffold,x=Mean),size=5,pch=21,fill="#98d5fa") +
  #scale_x_continuous(limits=c(5,35)) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.grid=element_blank()) + 
  scale_x_continuous(breaks=seq(5,75,5))

png(paste0(FIGURE_DIR,"/PAUL4_Depth_Autosomes_95&99Quaantiles.png"),res=300,width=7,height=6,units='in')
plot(p)
dev.off()

# This should give us enough information to filter on depth
