
# Plotting coverage statistics for WGS data

# Libraries

library(ggplot2)
library(scales)
library(tidyr)
library(dplyr)
library(patchwork)

# Set environment

source("/Volumes/Alter/Daus_WGS_Paper/config.R")
source(paste0(DATA_DIR,"/colours.R"))
setwd(paste0(WORKING_DIR,"/AnalysingWGSCoverage"))

# Get sample details
details <- read.table(paste0(HOME_DIR,"/References/wgs_sample_details.txt"),
                      header=T,stringsAsFactors=F,sep="\t")
details$Index <- 1:nrow(details)

# Rename populations
details$Pop <- factor(details$Pop,
                      levels=c("Wild","LHIP","LHISI"),
                      labels=c("Wild","Hybrid","Captive"))

# Get the file of colour schemes
colors <- read.delim(paste0(HOME_DIR,"/References/three_pop_colour_map.txt"),
                     header=T,stringsAsFactors=F,sep="\t")
cols <- colors[c(1,2,3),"Set1"]
names(cols) <- c("Wild","Hybrid","Captive")
my_colors <- scale_fill_manual(name = "Population", values = cols)
my_colors_col <- scale_color_manual(name = "Population", values = cols)

# Get actual depth data
data <- read.table("coverage_counts.txt", header=F,stringsAsFactors = F,sep="\t")
data <- data[grep("CM",data$V1),]
colnames(data) <- c("Chrom","Start","End","Mean","Sample")

# Merge
data <- merge(details,data)

# Plot mean + sd coverage

p1 <- data %>% group_by(Sample) %>%
  dplyr::summarise(Mean_Coverage = mean(Mean),
                   SD_Coverage = sd(Mean),
                   Sample=first(Sample),
                   Pop=first(Pop),
                   Sex=first(Sex)) %>%
  ggplot() + 
  geom_errorbar(aes(y=Mean_Coverage,
                    ymin=ifelse(Mean_Coverage - SD_Coverage < 0, 0, Mean_Coverage - SD_Coverage),
                    ymax=Mean_Coverage+SD_Coverage,
                    x=Sample),
                width=0) +
  geom_point(aes(y=Mean_Coverage,
                 x=Sample,
                 fill=Pop,
                 shape=Sex),size=3) + 
  scale_shape_manual(values=c(21,22)) +
  facet_grid(.~Pop,scales="free_x",space="free_x") + 
  theme_bw() +
  theme(axis.ticks=element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x=element_text(angle=30,hjust=1,vjust=1)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),limits=c(0,10),
                     breaks=c(0:10)) + 
  ylab("Mean depth of\ncoverage +/- SD") +
  my_colors + 
  guides(fill="none",shape="none")
  
p2 <- data %>% group_by(Sample) %>% summarise(mean_cov=mean(Mean),Sex=first(Sex),Pop=first(Pop)) %>%
  ggplot(aes(y=mean_cov)) + 
  geom_boxplot(outlier.alpha=0,width=0.05) + 
  geom_point(aes(x=0,fill=Pop,shape=Sex),position=position_jitter(width=0.01),size=3) + 
  theme_bw() +
  scale_shape_manual(values=c(21,22)) +
  theme(axis.ticks=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank()) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),limits=c(0,10),
                     breaks=c(0:10)) + 
  scale_x_continuous(expand=c(0.02,0.02)) +
  my_colors + 
  guides(fill="none")

p <- p1 + p2 + plot_layout(widths=c(8,1))

png(paste0(FIGURE_DIR,"/lcwgs_coverage_",format(Sys.time(),"%Y%m%d"),".png"),
    width=8,height=4,units='in',res=300)
plot(p)
dev.off()

# Now do the PAUL4 WGS sample

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

png(paste0(FIGURE_DIR,"/PAUL4_Depth_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=5,height=3,units='in')
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

png(paste0(FIGURE_DIR,"/PAUL4_Depth_Autosomes_95&99Quaantiles_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=7,height=6,units='in')
plot(p)
dev.off()

# This should give us enough information to filter on depth

# Also print out values for table

write.table(data %>%
  mutate(Type=ifelse(Chrom=="CM057009.1","Sex","Autosome")) %>%
  group_by(Sample,Type) %>%
  summarise(Mean=weighted.mean(Mean,End-Start)) %>%
  pivot_wider(names_from="Type",values_from="Mean"),
  "lcwgs_mean_coverage.txt",col.names=T,row.names=F,sep="\t",quote=F)

write.table(depths %>% 
  group_by(Type) %>%
  summarise(Mean=weighted.mean(Depth,End-Start)),
  "paul4_mean_coverage.txt",col.names=T,row.names=F,sep="\t",quote=F)

