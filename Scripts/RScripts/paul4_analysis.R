
### Script to analyse and plot depth in PAUL4 WGS sample

# Set environment

source("/Volumes/Alter/Daus_WGS_Paper/config.R")
setwd(paste0(WORKING_DIR,"/PAUL4"))

# Libraries

library(dplyr)
library(magrittr)
library(ggplot2)
library(Hmisc)
library(patchwork)
library(cowplot)
library(readxl,       warn.conflicts = F)
library(scales)
library(tidyr)
library(ggridges)

###### Coverage

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
  geom_histogram(colour="black",bins=40) + 
  scale_x_continuous(limits=c(5,30),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,7500),expand=c(0,0)) +
  scale_fill_manual(values=c("#98d5fa","#2e6094")) +
  theme_bw() +
  theme(axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.grid=element_blank())

png(paste0(FIGURE_DIR,"/PAUL4_Depth_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=6,height=3.5,units='in')
plot(p)
dev.off()

# Now analyse only autosomes
# How uniform is coverage? Can we use a universal filter?

grand_mean <- depths %>% 
  filter(Type == "Autosome") %>%
  summarise(Mean=mean(Depth)) %>%
  pull(Mean)

autosomes <- depths %>% 
  filter(Type == "Autosome") %>%
  mutate(Scaffold=as.factor(Scaffold))

p <- autosomes %>%
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
                width=0,linewidth=1) + 
  geom_errorbar(aes(y=Scaffold,x=Mean,xmin=Lower_99,xmax=Upper_99),
                width=0,linewidth=0.5) + 
  geom_point(aes(y=Scaffold,x=Mean),size=5,pch=21,fill="#98d5fa") +
  #scale_x_continuous(limits=c(5,35)) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.grid=element_blank()) + 
  scale_x_continuous(breaks=seq(5,75,5)) + 
  scale_y_discrete(limits=rev(levels(autosomes$Scaffold)),
                   labels=16:1,name="Chromosome"); p

png(paste0(FIGURE_DIR,"/PAUL4_Depth_Autosomes_95&99Quaantiles_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=7,height=6,units='in')
plot(p)
dev.off()

###### Genome wide H

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
colours <- rep(c("cornflowerblue","deepskyblue1"),8)

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
  summarise(W_Mean = weighted.mean(H,Coverage*(End-Start)),
            W_SE = weighted.se.mean(H,Coverage*(End-Start)))

# Make a final label to get genome wide average

genome_wide_h <- H %>%
  filter(Coverage > 0.1) %>%
  ggplot(aes(x=Bin,group=Chromosome,y=H,colour=Chromosome,fill=Chromosome)) +
  theme_bw() +
  theme(axis.ticks=element_blank(),panel.grid=element_blank(),
        legend.position="none",
        axis.text.x=element_text(hjust=1),
        axis.title.x=element_blank()) +
  geom_rect(inherit.aes=NULL,data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            colour="gray95",fill="gray95") +
  geom_line() +
  geom_area() +
  scale_x_continuous(expand=c(0,0),lim=c(1,nrow(H)),breaks=Locations$Locations,labels=1:16) +
  scale_y_continuous(breaks=y_labels,
                     labels=c(expression(0),
                              expression(2%*%10^-3),
                              expression(4%*%10^-3),
                              expression(6%*%10^-3),
                              expression(8%*%10^-3),
                              expression(1%*%10^-2))) + 
  coord_cartesian(ylim=c(ymin,ymax),expand=F) + 
  scale_colour_manual(unique(H$Chromosome),values=colours) +
  scale_fill_manual(unique(H$Chromosome),values=colours) +
  labs(y="Heterozygosity",x="Chromosome") + 
  geom_hline(yintercept=s_stats$W_Mean,linetype="dashed") ; genome_wide_h

png(paste0(FIGURE_DIR,"/genome_wide_H_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=10,height=5,units='in')
plot(genome_wide_h)
dev.off()

# Now we're going to look at the overall H for each chromosome
# We'll simply calculate a mean and a standard deviation, and plot this against chromosome length
# We take every second row to get non-overlapping windows

tmp <- H %>%
  filter(Coverage > 0.1) %>%
  group_by(Chromosome) %>%
  filter(row_number() %% 2 == 0) %>%
#  sample_n(100000,replace=T) %>%
  summarise(W_Mean = weighted.mean(H,Coverage),
            LQuantile = quantile(H,0.025),
            UQuantile = quantile(H,0.975),
            Length = first(max(End)))

tmp_model <- lm(W_Mean~Length,data=tmp)
summary(tmp_model)

H_v_length <- tmp %>%
  ggplot(aes(x=Length,y=W_Mean)) + 
  geom_smooth(method="lm",colour="black") +
  geom_point(size=3,pch=21,
             #aes(fill=Chromosome),
             fill="cornflowerblue",
             colour="black") + 
  #scale_fill_manual(unique(H$Chromosome),values=colours) +
  scale_y_continuous(limits=c(1.95e-4-0.00003,6.05e-4+0.00003),expand=c(0,0),
                     breaks=seq(from=2e-4,to=6e-4,by=1e-4),
                     labels=c(expression(2%*%10^-4),
                       expression(3%*%10^-4),
                       expression(4%*%10^-4),
                       expression(5%*%10^-4),
                       expression(6%*%10^-4))) + 
  scale_x_continuous(breaks=c(1e8,2e8,3e8,4e8),
                     labels=c(1,2,3,4)) +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.grid=element_blank(),
        legend.position="none") + 
  labs(x=expression(paste("Chromosome length (",bp%*%10^8,")")),
       y="Heterozygosity") + 
  annotate(geom="text",
           x = Inf, y = Inf,
           hjust = 1.5, vjust = 4.5,
           label=paste0("p = ",
                        round(summary(tmp_model)$coefficients[2,4],3))); H_v_length

png(paste0(FIGURE_DIR,"/H_v_chrom_length_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=5,height=5,units='in')
plot(H_v_length)
dev.off()

# Now we're going to add to this the figure of genome_wide H of island species
# This is from published data on the Kakapo and the island foxes

islands <- suppressWarnings(read_excel(paste0(DATA_DIR,"/20231222_published_h_estimates.xlsx"),sheet="Sheet1"))

islands <- rbind(islands,
      c("Dryococelus_australis","LHISI","D. australis",s_stats$W_Mean,"This study",NA)) %>%
  filter(Clade %in% c("Chordata","Insecta","D. australis")) %>%
  mutate(Observed_H = as.numeric(Observed_H)) %>%
  group_by(Common_Name) %>%
  summarise(Observed_H = mean(Observed_H),Clade=first(Clade)) %>%
  ungroup()

arrow_position=s_stats$W_Mean-0.000025
published_h <- islands %>%
  ggplot(aes(x=Observed_H,fill=Clade)) + 
  geom_histogram(colour="black",bins=30) + 
  scale_y_continuous(limits=c(0,6.5),expand=c(0,0)) +
  scale_fill_manual(breaks=c("Chordata","Insecta","D. australis"),
                    labels=c("Chordata","Insecta",expression(italic("D. australis"))),
                    values=c("grey88","grey64","cornflowerblue")) + 
  scale_x_continuous(trans="log10",
                     breaks=c(2e-4,2e-3,2e-2),
                     labels=c(expression(2%*%10^-4),
                              expression(2%*%10^-3),
                              expression(2%*%10^-2))) +
  # scale_x_continuous(trans="log10",
  #                    breaks=c(1e-4,1e-3,1e-2),
  #                    labels=c(expression(1%*%10^-4),
  #                             expression(1%*%10^-3),
  #                             expression(1%*%10^-2))) +
  theme_bw() + 
  theme(
#        axis.text.y=element_blank(),
#        axis.title.y=element_blank(),
        axis.ticks=element_blank(),
        legend.position=c(0.85,0.85),
        legend.title = element_blank(),
        legend.background = element_rect(colour="black"),
        panel.grid=element_blank()) + 
  labs(x="Heterozygosity",y="Count") + 
  annotate(geom="segment",
           y=3,yend=2.2,x=arrow_position,xend=arrow_position,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) ; published_h

png(paste0(FIGURE_DIR,"/published_genomic_H_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=5,height=5,units='in')
plot(published_h)
dev.off()

# Also write list of species with lower H

islands %>% filter(Observed_H <= s_stats$W_Mean)
# Baiji, Cheetah, Island_fox, Snow_leopard, Tasmanian_devil

# Finally, make a plot of per chromosome GC content

paul4 <- read.table("paul4_gc.txt",header=T,stringsAsFactors = F)
lengths <- read.table(paste0(REF_DIR,"/scaffold_lengths"),header=F,stringsAsFactors = F)
colnames(lengths) <- c("chr","total_length")

paul4_gc <- paul4 %>% 
  group_by(chr) %>% 
#  sample_n(10000,replace=T) %>% 
  mutate(GC = (G+C)/length) %>%
  dplyr::summarise(mean_GC=mean(GC)) %>%
  merge(.,lengths) 

tmp_model <- lm(mean_GC~total_length,data=paul4_gc)
summary(tmp_model)

paul4_gc_plot <- paul4_gc %>%
  ggplot(aes(x=total_length,y=mean_GC*100)) + 
  geom_smooth(method="lm",colour="black") +
  geom_point(size=3,pch=21,
             #aes(fill=chr),
             fill="cornflowerblue",
             colour="black") + 
  #scale_fill_manual(unique(paul4_gc$chr),values=colours) +
  scale_x_continuous(breaks=c(1e8,2e8,3e8,4e8),
                     labels=c(1,2,3,4)) +
  scale_y_continuous(limits=c(38.1,39.25),expand=c(0,0),
                     breaks=seq(from=38.2,to=39.2,0.2)) +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.grid=element_blank(),
        legend.position="none") + 
  labs(x=expression(paste("Chromosome length (",bp%*%10^8,")")),
       y="GC content (%)") +
  annotate(geom="text",
           x = Inf, y = Inf,
           hjust = 1.5, vjust = 4.5,
           label=paste0("p = ",
                        round(summary(tmp_model)$coefficients[2,4],3)))

png(paste0(FIGURE_DIR,"/paul_gc_content_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=5,height=5,units='in')
plot(paul4_gc_plot)
dev.off()

# Also save all of these as one figure

layout <- "
AAAA
AAAA
BBDD
CCDD
"
paul4_summary <- genome_wide_h + 
  (H_v_length + theme(axis.text.x=element_blank())) + 
  paul4_gc_plot + 
  published_h + 
  plot_layout(design = layout,axis_titles = "collect")

png(paste0(FIGURE_DIR,"/paul_summary_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=10.5,height=10.5,units='in')
plot(paul4_summary)
dev.off()

# Is there a relationship between H and the distance from the centre of the chromosome?
p1 <- H %>% 
  group_by(Chromosome) %>%
  mutate(Midpoint=(End+Start)/2,
         Distance=abs(Midpoint-median(Midpoint))/max(Midpoint)) %>%
  ggplot(aes(x=Distance,y=log(H))) + 
  geom_point() + 
  geom_smooth(method="lm") + 
  facet_wrap(~Chromosome) +
  theme_bw() + 
  theme(axis.ticks=element_blank())
p2 <- H %>% 
  group_by(Chromosome) %>%
  mutate(Midpoint=(End+Start)/2,
         Distance=abs(Midpoint-median(Midpoint))/max(Midpoint)) %>%
  ggplot(aes(x=Distance,y=log(H))) + 
  geom_point() + 
  geom_smooth(method="lm") +
  theme_bw() + 
  theme(axis.ticks=element_blank())
p1 + p2 + plot_layout(widths=c(6,4))
# In general, yes

# One final figure...

gc_bic <- read.table("paul4_gc_big_windows.txt",header=T)
colnames(gc_bic) <- c("Chromosome","Start","End","A","C","G","T")
gc_bic <- merge(gc_bic,H)
gc_bic <- gc_bic %>%
  mutate(GC = (G+C)/(End-Start))

gc_bic %>%
  ggplot(aes(x=GC,y=H)) + 
  geom_point() + 
  facet_wrap(~Chromosome) + 
  geom_smooth(method="lm")

##### Genome wide ROH

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
library(readxl,       warn.conflicts = F)

# Get lengths
lengths <- read.table("../../References/scaffold_lengths")[1:16,]
colnames(lengths) <- c("CHR","Total")

# BCFTOOLS

roh_bcf_04 <- readr::read_table("GENOME_WIDE_ROH_04_BCFTOOLS.txt") ; colnames(roh_bcf_04)[c(1,4)] <- c("CHR","Length")
roh_bcf_05 <- readr::read_table("GENOME_WIDE_ROH_05_BCFTOOLS.txt") ; colnames(roh_bcf_05)[c(1,4)] <- c("CHR","Length")
roh_bcf_06 <- readr::read_table("GENOME_WIDE_ROH_06_BCFTOOLS.txt") ; colnames(roh_bcf_06)[c(1,4)] <- c("CHR","Length")

### Get a single data.frame with all lengths

roh_bcf_04$AF = "AF04" ; roh_bcf_05$AF = "AF05" ; roh_bcf_06$AF = "AF06"
all_roh <- rbind(roh_bcf_04,roh_bcf_05,roh_bcf_06) %>% merge(.,lengths)

# Which chromosome has the greatest proportion of homozygosity
all_roh %>%
  filter(AF=="AF05") %>%
  group_by(CHR) %>%
  summarise(FROH=sum(Length)/first(Total)) %>%
  arrange(-FROH)
# Which chromosome had the longest
all_roh %>%
  filter(AF=="AF05") %>%
  group_by(CHR) %>%
  filter(Length==max(Length)) %>%
  ungroup() %>%
  arrange(-Length)
# Do chromosomes differ in their average lengths
all_roh %>%
  filter(AF=="AF05") %>%
  group_by(CHR) %>%
  summarise(Mean=mean(Length)) %>%
  arrange(Mean)
all_roh %>%
  filter(AF=="AF05") %>%
  group_by(CHR) %>%
  summarise(Mean=mean(Length)) %>%
  ungroup() %>%
  summarise(Mean_Mean=mean(Mean),SD=sd(Mean))
# Where are the longest ROH on the two chromosomes with the longest ROH
all_roh %>%
  filter(AF=="AF05",CHR %in% c("CM056996.1","CM056995.1")) %>%
  group_by(CHR) %>%
  arrange(-Length) %>%
  slice_head(n=5)

# Calculate ROH using different length cutoffs

roh_filt <- function(cutoff = 0,dataset=NULL){
  
  if(is.null(data)){stop("supply a data.frame")}
  
  dataset %>%
    filter(Length > cutoff) %>%
    group_by(CHR) %>%
    summarise(Length = sum(Length)) %>%
    summarise(sum(Length)/sum(lengths$Total)) %>%
    return
  
}

cutoffs <- c(1e5,1e6,2e6,5e6)

tmp <- data.frame(ROH=c(unlist(lapply(cutoffs,FUN=roh_filt,dataset=roh_bcf_04)),
                 unlist(lapply(cutoffs,FUN=roh_filt,dataset=roh_bcf_05)),
                 unlist(lapply(cutoffs,FUN=roh_filt,dataset=roh_bcf_06))),
           Length = rep(cutoffs,3),
           Method = rep(c("AF = 0.4",
                          "AF = 0.5",
                          "AF = 0.6"),each=4))

ROHs_figure_1 <- tmp %>%
  #filter(Method != "PLINK") %>%
  ggplot(aes(x=as.factor(Length),y=ROH,fill=Method)) + 
  geom_bar(stat="identity",position=position_dodge(),colour="black") + 
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.grid=element_blank(),
        legend.position.inside=c(0.85,0.86),
        legend.background = element_rect(colour="black"),
        legend.title = element_blank()) + 
  scale_x_discrete(labels=c("> 100 kb","> 1 Mb","> 2 Mb","> 5 Mb")) + 
  labs(x="Length cutoff",y=expression(F[ROH])) + 
  scale_fill_brewer(palette = "Blues")

png(paste0(FIGURE_DIR,"/roh_comparison_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=5.5,height=5.5,units='in')
plot(ROHs_figure_1)
dev.off()

# Also save this as a table

write.table(tmp,paste0("roh_comparison_",format(Sys.time(),"%Y%m%d"),".txt"),
            row.names=F,col.names=T,quote=F,sep="\t")

### Count of ROH per size category

count_function <- function(data=NULL,label_=NULL){
  
  if(is.null(data)){stop("supply a data.frame")}
  
  data$Bin <- cut(data$Length,c(1e5,1e6,2e6,3e6,4e6,5e6,100e6))
  data %>% group_by(Bin) %>%
    summarise(Count = n()) %>%
    filter(!is.na(Bin)) %>%
    mutate(Method=label_) %>%
    return
  
}

tmp <- count_function(roh_bcf_04,label_="AF = 0.4") %>%
  merge(.,count_function(roh_bcf_05,label_="AF = 0.5"),all=T) %>%
  merge(.,count_function(roh_bcf_06,label_="AF = 0.6"),all=T) %>% 
  tidyr::complete(.,Bin,Method) %>%
  tidyr::replace_na(list(Count=0))

tmp

ROHs_figure_2 <- tmp %>%
  #filter(Method != "PLINK") %>%
  ggplot(aes(x=Bin,y=Count,fill=Method)) + 
  geom_bar(stat="identity",position=position_dodge(),colour="black") + 
  scale_y_continuous(limits=c(0,820),expand=c(0,0)) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.grid=element_blank(),
        legend.position=c(0.85,0.86),
        legend.background = element_rect(colour="black"),
        legend.title=element_blank()) + 
  scale_x_discrete(labels=c("100 kb -\n1 Mb","1 Mb -\n2 Mb",
                            "2 Mb -\n3 Mb","3 Mb -\n4 Mb","4 Mb -\n5 Mb",
                            "> 5 Mb")) + 
  labs(x="Length",y="Count") + 
  scale_fill_brewer(palette = "Blues")

png(paste0(FIGURE_DIR,"/roh_count_comparison_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=5.5,height=5.5,units='in')
plot(ROHs_figure_2)
dev.off()

# Also save as table

write.table(tmp,paste0("roh_count_comparison_",format(Sys.time(),"%Y%m%d"),".txt"),
            row.names=F,col.names=T,quote=F,sep="\t")

### Get table of stats for all methods
### Per chrom + overall mean, median, min, max, sd, etc

tmp <- roh_plink %>%
  filter(Length > 1e5) %>% group_by(CHR) %>%
  summarise(Mean=mean(Length),
            SD=sd(Length),
            Max=max(Length),
            Min=min(Length),
            Median=median(Length),
            Count=n(),
            Method="PLINK")
tmp1 <- roh_plink %>% filter(Length > 1e5) %>%
  summarise(CHR="Overall",
            Mean=mean(Length),
            SD=sd(Length),
            Max=max(Length),
            Min=min(Length),
            Median=median(Length),
            Count=n(),
            Method="PLINK")
roh_plink_stats <- rbind(tmp,tmp1)

tmp <- roh_bcf_04 %>%
  filter(Length > 1e5) %>% group_by(CHR) %>%
  summarise(Mean=mean(Length),
            SD=sd(Length),
            Max=max(Length),
            Min=min(Length),
            Median=median(Length),
            Count=n(),
            Method="BCF_AF_04")
tmp1 <- roh_bcf_04 %>% filter(Length > 1e5) %>%
  summarise(CHR="Overall",
            Mean=mean(Length),
            SD=sd(Length),
            Max=max(Length),
            Min=min(Length),
            Median=median(Length),
            Count=n(),
            Method="BCF_AF_04")
roh_bcf_04_stats <- rbind(tmp,tmp1)

tmp <- roh_bcf_05 %>%
  filter(Length > 1e5) %>% group_by(CHR) %>%
  summarise(Mean=mean(Length),
            SD=sd(Length),
            Max=max(Length),
            Min=min(Length),
            Median=median(Length),
            Count=n(),
            Method="BCF_AF_05")
tmp1 <- roh_bcf_05 %>% filter(Length > 1e5) %>%
  summarise(CHR="Overall",
            Mean=mean(Length),
            SD=sd(Length),
            Max=max(Length),
            Min=min(Length),
            Median=median(Length),
            Count=n(),
            Method="BCF_AF_05")
roh_bcf_05_stats <- rbind(tmp,tmp1)

tmp <- roh_bcf_06 %>%
  filter(Length > 1e5) %>% group_by(CHR) %>%
  summarise(Mean=mean(Length),
            SD=sd(Length),
            Max=max(Length),
            Min=min(Length),
            Median=median(Length),
            Count=n(),
            Method="BCF_AF_06")
tmp1 <- roh_bcf_06 %>% filter(Length > 1e5) %>%
  summarise(CHR="Overall",
            Mean=mean(Length),
            SD=sd(Length),
            Max=max(Length),
            Min=min(Length),
            Median=median(Length),
            Count=n(),
            Method="BCF_AF_06")
roh_bcf_06_stats <- rbind(tmp,tmp1)

write.table(roh_plink_stats,paste0("roh_plink_stats_",format(Sys.time(),"%Y%m%d"),".txt"),
            row.names=F,col.names=T,quote=F,sep="\t")
write.table(roh_bcf_04_stats,paste0("roh_bcf_04_stats_",format(Sys.time(),"%Y%m%d"),".txt"),
            row.names=F,col.names=T,quote=F,sep="\t")
write.table(roh_bcf_05_stats,paste0("roh_bcf_05_stats_",format(Sys.time(),"%Y%m%d"),".txt"),
            row.names=F,col.names=T,quote=F,sep="\t")
write.table(roh_bcf_06_stats,paste0("roh_bcf_06_stats_",format(Sys.time(),"%Y%m%d"),".txt"),
            row.names=F,col.names=T,quote=F,sep="\t")

# Plotting ROH length distribution in detail

breaks <- c(c(1,2,5) * rep(c(1e4,1e5,1e6),each=3),1e7)

p <- all_roh %>%
  ggplot(aes(x=Length)) +
  geom_histogram(aes(fill=AF),colour="black") + 
  scale_fill_brewer(palette="Blues",guide="none") + 
  scale_x_continuous(trans="log10",
                     breaks=breaks) +
  facet_wrap(~AF,ncol=1) + 
  theme_bw() + 
  ylab("Count") ; p

png(paste0(FIGURE_DIR,"/roh_lengths_",format(Sys.time(),"%Y%m%d"),".png"),
    res=300,width=5.5,height=5.5,units='in')
plot(p)
dev.off()

all_roh %>%
  group_by(AF) %>%
  dplyr::summarise(Mean=mean(Length),
                   SE=sd(Length)/sqrt(n()),
                   Max=max(Length))
