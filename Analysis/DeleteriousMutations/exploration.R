# Set directory structure
HOME_DIR="/Volumes/Alter/LHISI"
REF_DIR=paste0(HOME_DIR,"/References")
WORKING_DIR=paste0(HOME_DIR,"/Analyses/DeleteriousMutations")

setwd(WORKING_DIR)

# Load libraries
library(ggplot2)
library(dplyr)
library(plyr)
library(ggvenn)
library(tidyr)
library(ggxmean)

# Function for complement
`%ni%` <- Negate(`%in%`)

### Loading and diagnostics

# Read data in, rename for ease of use
lhisi <- read.table("../AlleleFrequencies/lhisi.mafs.gz",header=T,stringsAsFactors=F,sep="\t")[,c(1,2,6,7)]
lhip <- read.table("../AlleleFrequencies/lhip.mafs.gz",header=T,stringsAsFactors=F,sep="\t")[,c(1,2,6,7)]
wild <- read.table("../AlleleFrequencies/wild.mafs.gz",header=T,stringsAsFactors=F,sep="\t")[,c(1,2,6,7)]
colnames(lhisi) <- c("Scaffold","Position","MAF_lhisi","N_lhisi")
colnames(lhip) <- c("Scaffold","Position","MAF_lhip","N_lhip")
colnames(wild) <- c("Scaffold","Position","MAF_wild","N_wild")

# Read in effect data
effects <- read.table("vars_classified.txt",
                      header=T,stringsAsFactors=F,sep="\t")[,c(1,2,3,4,5,6,7,8,10)]

# Get site information
scafs <- c(1:17)[-4]
sites <- data.frame(Scaffold=character(),Position=character(),Ref=character(),Alt=character())
for(i in 1:length(scafs)){
  tmp <- eval(parse(text=paste0("read.table('../AlleleFrequencies/scaffold",scafs[i],"_sites',header=F,stringsAsFactors=F)")))
  sites <- rbind(sites,tmp)
}
rm(tmp) ; colnames(sites) <- c("Scaffold","Position","Ref","Alt")

# Give all sites their alleles
effects <- merge(effects,sites)
rm(sites)

# Get gene AED scores
scores <- read.table("gene_scores.txt",header=T,stringsAsFactors=F,sep="\t")[,c(1,2)]

# Give variants inside genes their AED score
intergenic <- effects %>% filter(is.na(Gene) == T)
genic <- effects %>% filter(is.na(Gene) == F)
genic <- merge(x=genic,y=scores,all.x=T)
effects <- rbind.fill(genic,intergenic)
rm(genic,intergenic,scores)

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

effects <- effects %>% mutate(captive_fate=ifelse((fate_wild == "absent" & 
                                                     fate_lhisi != "absent" &  
                                                     fate_lhip != "absent"), "common",
                                       ifelse((fate_wild == "absent" &  
                                                 fate_lhisi != "absent" &  
                                                 fate_lhip == "absent"), "lhisi",
                                              ifelse((fate_wild == "absent" &  
                                                        fate_lhisi == "absent" &  
                                                        fate_lhip != "absent"),  "lhip",NA))))

### Classifying variants based on what happens to them in captivity
### Generating more summary variables

# Function to get categorical variable of sign

get_sign <- function(x){
  # If sign is positive, decrease in captivity
  if( x > 0 ){return("increase")} else
    # If sign is negative, increase in captivity
    if ( x < 0 ){return("decrease")} else
      # If sign is NA, no change, fixed allele remains fixed
      if ( x == 0 ){return("no change")} else
        # Else
        if( is.na(x) ){return(NA)}
}


wild <- effects %>% filter(fate_wild == "segregating")

wild$lhisi_change <- wild$MAF_lhisi - wild$MAF_wild
wild$lhip_change <- wild$MAF_lhip - wild$MAF_wild

wild <- wild %>% drop_na(lhisi_change,lhip_change)

wild$lhisi_sign <- sapply(wild$lhisi_change,get_sign)
wild$lhip_sign <- sapply(wild$lhip_change,get_sign)

# Classifying variants based on presence in putative ROH in the three populations

lhip <- read.table("LHIP.cov",header=F,stringsAsFactors=F,sep="\t")
lhisi <- read.table("LHISI.cov",header=F,stringsAsFactors=F,sep="\t")
wild <- read.table("WILD.cov",header=F,stringsAsFactors=F,sep="\t")
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

# Plots

effects %>% 
#  mutate(new_MAF_lhisi = ifelse(MAF_lhisi > 0.5, 1-MAF_lhisi, MAF_lhisi)) %>%
  filter(fate_lhisi != "absent",Variant_Effect != "MODIFIER") %>% 
  ggplot(aes(x=MAF_lhisi,y=Coverage_lhisi)) + 
  facet_wrap(~Variant_Effect) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour="black",fill=NA)) + 
  stat_density_2d(
    geom = "raster",
    aes(fill = after_stat(sqrt(density))),
    contour = FALSE) + 
  stat_density_2d(colour="black") +
  #  scale_fill_viridis_c() +
  #  scale_fill_gradient(low="yellow",high="black") +
  scale_fill_gradientn(colours = terrain.colors(9)) +
  scale_x_continuous(limits=c(0,1),expand=c(0,0)) + 
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) + 
  labs(x="Minor allele frequency",y="ROH frequency")

effects %>% filter(fate_lhip != "absent",Variant_Effect != "MODIFIER") %>% 
  ggplot(aes(x=MAF_lhip,y=Coverage_lhip)) + 
  facet_wrap(~Variant_Effect) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour="black",fill=NA)) + 
  stat_density_2d(
    geom = "raster",
    aes(fill = after_stat(sqrt(density))),
    contour = FALSE) + 
  stat_density_2d(colour="black") +
#  scale_fill_viridis_c() +
#  scale_fill_gradient(low="yellow",high="black") +
  scale_fill_gradientn(colours = terrain.colors(9)) +
  scale_x_continuous(limits=c(0,1),expand=c(0,0)) + 
  scale_y_continuous(limits=c(0,1),expand=c(0,0)) + 
  labs(x="Minor allele frequency",y="ROH frequency")

# Now make some plots?

wild <- effects %>% filter(fate_wild != "absent")

wild$lhisi_change <- wild$MAF_lhisi - wild$MAF_wild
wild$lhip_change <- wild$MAF_lhip - wild$MAF_wild

wild <- wild %>% drop_na(lhisi_change,lhip_change)

wild$lhisi_sign <- sapply(wild$lhisi_change,get_sign)
wild$lhip_sign <- sapply(wild$lhip_change,get_sign)

ggplot(wild,aes(x=lhisi_change,fill=fate_lhisi)) + 
  geom_histogram(position="identity",colour="black") + 
  facet_wrap(~Variant_Effect,scales="free_y",ncol=1) +
  geom_x_median(size=1.5) + 
  scale_y_continuous(limits=c(0,NA),expand=c(0,0))

ggplot(wild,aes(x=lhip_change,fill=fate_lhip)) + 
  geom_histogram(position="identity",colour="black") + 
  facet_wrap(~Variant_Effect,scales="free_y",ncol=1) +
  geom_x_median(size=1.5) + 
  scale_y_continuous(limits=c(0,NA),expand=c(0,0))

wild %>% group_by(Variant_Effect) %>% 
  dplyr::summarise(mean_lhip=mean(lhip_change),
                   sd_lhip=sd(lhip_change),
                   mean_lhisi=mean(lhisi_change),
                   sd_lhisi=sd(lhisi_change)) %>%
  gather(key, value, -Variant_Effect) %>%
  extract(key, c("stat", "pop"), "(.+)_(.+)") %>%
  spread(stat, value) %>% 
  ggplot(aes(x=Variant_Effect,y=mean,fill=pop)) + 
  geom_point(size=3,pch=21,position=position_dodge(width=0.3)) + 
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0,position=position_dodge(width=0.3)) + 
  scale_y_continuous(limits=c(-0.8,0.1)) + 
  geom_hline(yintercept=0,linetype="dashed",colour="black",size=1.5) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.minor=element_blank()) + 
  labs(x="Variant\neffect",y="Frequency difference\nin captivity")

p1 <- wild %>% group_by(Variant_Effect) %>% 
  dplyr::summarise(mean=mean(lhip_change),
                   sd=sd(lhip_change),
                   med=median(lhip_change)) %>% 
  ggplot(aes(x=Variant_Effect,y=mean)) + 
  geom_point(size=3,pch=21) + 
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0) + 
  scale_y_continuous(limits=c(-0.8,0.1)) + 
  geom_hline(yintercept=0,linetype="dashed",colour="black",size=1.5) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.minor=element_blank()) + 
  labs(x="Variant\neffect",y="Frequency difference\nin captivity")

p2 <- wild %>% group_by(Variant_Effect) %>% 
  dplyr::summarise(mean=mean(lhisi_change),
                   sd=sd(lhisi_change),
                   med=median(lhisi_change)) %>% 
  ggplot(aes(x=Variant_Effect,y=mean)) + 
  geom_point(size=3,pch=21) + 
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0) + 
  scale_y_continuous(limits=c(-0.8,0.1)) + 
  geom_hline(yintercept=0,linetype="dashed",colour="black",size=1.5) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.minor=element_blank()) + 
  labs(x="Variant\neffect",y="Frequency difference\nin captivity")

p1 + p2

# Different aggregation

wild <- effects %>% filter(fate_wild != "absent")
wild$lhisi_change <- wild$MAF_lhisi - wild$MAF_wild
wild$lhip_change <- wild$MAF_lhip - wild$MAF_wild
wild <- wild %>% drop_na(lhisi_change,lhip_change)
wild$lhisi_sign <- sapply(wild$lhisi_change,get_sign)
wild$lhip_sign <- sapply(wild$lhip_change,get_sign)

wild %>% ggplot(aes(x=lhisi_change,fill=fate_wild)) + 
  geom_histogram(position="identity",colour="black") + 
  facet_wrap(~Variant_Effect,scales="free_y",ncol=1) +
  geom_x_median(size=1.5) + 
  scale_y_continuous(limits=c(0,NA),expand=c(0,0))

wild %>% ggplot(aes(x=lhip_change,fill=fate_wild)) + 
  geom_histogram(position="identity",colour="black") + 
  facet_wrap(~Variant_Effect,scales="free_y",ncol=1) +
  geom_x_median(size=1.5) + 
  scale_y_continuous(limits=c(0,NA),expand=c(0,0))

wild %>% ggplot(aes(x=lhisi_change,fill=fate_lhisi)) + 
  geom_histogram(position="identity",colour="black") + 
  facet_wrap(~Variant_Effect+fate_wild,scales="free_y",ncol=2) +
  geom_x_median(size=1.5) + 
  scale_y_continuous(limits=c(0,NA),expand=c(0,0)) + 
  scale_x_continuous(limits=c(-1.05,1.05)) + 
  geom_vline(xintercept=0,linetype="dashed")

wild %>% ggplot(aes(x=lhip_change,fill=fate_lhip)) + 
  geom_histogram(position="identity",colour="black") + 
  facet_wrap(~Variant_Effect+fate_wild,scales="free_y",ncol=2) +
  geom_x_median(size=1.5) + 
  scale_y_continuous(limits=c(0,NA),expand=c(0,0)) + 
  scale_x_continuous(limits=c(-1.05,1.05)) + 
  geom_vline(xintercept=0,linetype="dashed")


