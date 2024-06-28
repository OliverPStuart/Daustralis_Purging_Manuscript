### Script to analyse and plot H from lcWGS

# Libraries

library(ggplot2)
library(dplyr)
library(plyr)
library(tidyr)
library(ggnewscale)
library(ggResidpanel)
options(dplyr.summarise.inform = FALSE)

# Set environment

source("../../config.R")
setwd(paste0(WORKING_DIR,"/DeleteriousMutations"))

# My favourite function

`%ni%` <- Negate(`%in%`)

# Get colours

source(paste0(DATA_DIR,"/colours.R"))

# Read data in, rename for ease of use

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
effects <- rbind.fill(genic,intergenic)
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

# Bring in 2-dimensional colour scheme and make legend figures

new_colours <- read.delim(paste0(REF_DIR,"/two_pop_2d_colour_scale.txt"),
                          header=F,stringsAsFactors = F)
colnames(new_colours) <- c("Effect","Captive","Hybrid","Wild")

new_colours_n <- new_colours %>% gather(-Effect,key="Population",value="Colour")

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

# Get wild alleles

wild <- effects %>% filter(fate_wild != "absent" & !is.na(fate_wild))

# Calculate difference in frequency

wild$lhisi_change <- wild$MAF_lhisi - wild$MAF_wild
wild$lhip_change <- wild$MAF_lhip - wild$MAF_wild
wild <- wild %>% drop_na(lhisi_change,lhip_change)

# Plotting mean+sd difference per variant effect class
# Excluding those variants in each population which are absent

tmp_lhisi <- wild %>% group_by(Variant_Effect) %>% 
#  filter(is.na(AED) | AED > 0.25) %>%
  filter(fate_lhisi != "absent") %>%
  dplyr::summarise(mean_lhisi=mean(lhisi_change),
                   sd_lhisi=sd(lhisi_change)) %>%
  gather(key, value, -Variant_Effect) %>%
  tidyr::extract(key, c("stat", "pop"), "(.+)_(.+)") %>%
  spread(stat, value)
tmp_lhip <- wild %>% group_by(Variant_Effect) %>% 
  #  filter(fate_lhip != "absent") %>%
#  filter(is.na(AED) | AED > 0.25) %>%
  dplyr::summarise(mean_lhip=mean(lhip_change),
                   sd_lhip=sd(lhip_change)) %>%
  gather(key, value, -Variant_Effect) %>%
  tidyr::extract(key, c("stat", "pop"), "(.+)_(.+)") %>%
  spread(stat, value)

change_data <- rbind(tmp_lhisi,tmp_lhip)

p <- change_data %>% 
  ggplot() + 
  scale_x_continuous(limits=c(-0.7,0.1)) + 
  geom_errorbar(aes(fill=factor(Variant_Effect,levels=c("HIGH","MODERATE","LOW","MODIFIER")),
                    x=mean,y=pop,xmin=mean-sd,xmax=mean+sd),
                width=0,size=0.7,position=position_dodge(width=0.6),
                filter(change_data,pop=="lhisi")) + 
  geom_point(aes(fill=factor(Variant_Effect,levels=c("HIGH","MODERATE","LOW","MODIFIER")),
                 x=mean,y=pop),
             size=4.5,pch=21,position=position_dodge(width=0.6),
             filter(change_data,pop=="lhisi")) + 
  scale_fill_manual(values=rev(new_colours[,3])) +
  new_scale_fill() +
  geom_errorbar(aes(fill=factor(Variant_Effect,levels=c("HIGH","MODERATE","LOW","MODIFIER")),
                    x=mean,y=pop,xmin=mean-sd,xmax=mean+sd),
                width=0,size=0.7,position=position_dodge(width=0.6),
                filter(change_data,pop=="lhip")) + 
  geom_point(aes(fill=factor(Variant_Effect,levels=c("HIGH","MODERATE","LOW","MODIFIER")),
                 x=mean,y=pop),
             size=4.5,pch=21,position=position_dodge(width=0.6),
             filter(change_data,pop=="lhip")) + 
  scale_fill_manual(values=rev(new_colours[,2])) +
  geom_vline(xintercept=0,linetype="dashed",colour="black",size=1) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        panel.grid=element_blank(),
        legend.position="none") + 
  labs(x="Frequency difference\nin captivity") +
  scale_y_discrete(labels=c("LHIP","LHISI")) ; p

# Plot

png(paste0(FIGURE_DIR,"/frequency_change_",format(Sys.time(),"%Y%m%d"),".png"),
    width=6,height=3,units='in',res=300)
plot(p)
dev.off()

# Now run the tests
# We exclude WildPresent-CaptiveAbsent alleles here
# However, we include WildPresent-HybridAbsent

tmp_lhisi <- wild %>% filter(fate_lhisi != "absent") %>% 
#  filter(is.na(AED) | AED > 0.25) %>%
  dplyr::select(Variant_Effect,lhisi_change) %>% mutate(Pop="LHISI")
tmp_lhip <- wild %>%
#  filter(is.na(AED) | AED > 0.25) %>%
  dplyr::select(Variant_Effect,lhip_change) %>% mutate(Pop="LHIP")
colnames(tmp_lhip)[2] <- "change" ; colnames(tmp_lhisi)[2] <- "change"
tmp <- rbind(tmp_lhisi,tmp_lhip)

lhisi_model <- lm(change~Variant_Effect,data=tmp[tmp$Pop=="LHISI",])
lhip_model <- lm(change~Variant_Effect,data=tmp[tmp$Pop=="LHIP",])

#resid_panel(lhisi_model)
#resid_panel(lhip_model)

summary(lhisi_model)
summary(lhip_model)

# Ratio of deleterious variant presence inside/outside of high ROH regions
# Need to standardise this by presence of ROH in the first place....
# Would that just alter the difference between populations and not the populations themselves?

# Getting some sampling error intervals for these with bootstrapping

lhisi_base <- effects %>% filter(
  fate_lhisi == "segregating",
  #Variant_Effect != "MODIFIER"
) %>%
#  filter(is.na(AED) | AED > 0.25) %>%
  mutate(ROH_50 = ifelse(Coverage_lhisi > 0.5,"Upper","Lower"))

lhip_base <- effects %>% filter(
  fate_lhip == "segregating"
) %>%
#  filter(is.na(AED) | AED > 0.25) %>%
  mutate(ROH_50 = ifelse(Coverage_lhip > 0.5,"Upper","Lower"))

# Now bootstrap

# Number of bootstraps
boot=100
# Get variant vector
types <- c("MODIFIER","LOW","MODERATE","HIGH")

# LHISI
vec <- c()
for(i in 1:boot){
  
  #if(i %% 100 == 0){
    print(i)
  #}
  
  val <- lhisi_base[sample(x=1:nrow(lhisi_base),size=nrow(lhisi_base),replace=T),] %>%
    group_by(ROH_50,Variant_Effect) %>%
    dplyr::summarise(tabulate=n()) %>% 
    spread(key="ROH_50",value="tabulate") %>% 
    mutate(ratio = Lower / (Lower+Upper)) %>%
    pull(ratio)
  
  vec <- c(vec,val)
  
  if(i == boot){
    bootstrap_lhisi <- data.frame(ratio=vec,
                                  Variant_Effect=rep(types,boot)) %>%
      group_by(Variant_Effect) %>%
      dplyr::summarise(sd=sd(ratio),
                       mean=mean(ratio),
                       median=median(ratio),
                       upper_conf_95=quantile(ratio, probs = c(0.025,0.975))[2],
                       lower_conf_95=quantile(ratio, probs = c(0.025,0.975))[1],
                       Pop="LHISI")
    
    bootstrap_lhisi$Variant_Effect <- factor(bootstrap_lhisi$Variant_Effect,
                                             labels=types,
                                             levels=types)
    
    ratio <- lhisi_base %>%
      group_by(Variant_Effect,ROH_50) %>%
      dplyr::summarise(tabulate=n()) %>% 
      spread(key="ROH_50",value="tabulate") %>% 
      mutate(ratio = Lower / (Lower+Upper))
    
    bootstrap_lhisi <- merge(ratio,bootstrap_lhisi)
    
  }
  
}

# LHIP
vec <- c()
for(i in 1:boot){
  
  #if(i %% 100 == 0){
    print(i)
  #}
  
  val <- lhip_base[sample(x=1:nrow(lhip_base),size=nrow(lhip_base),replace=T),] %>%
    group_by(ROH_50,Variant_Effect) %>%
    dplyr::summarise(tabulate=n()) %>% 
    spread(key="ROH_50",value="tabulate") %>% 
    mutate(ratio = Lower / (Lower+Upper)) %>%
    pull(ratio)
  
  vec <- c(vec,val)
  
  if(i == boot){
    bootstrap_lhip <- data.frame(ratio=vec,
                                 Variant_Effect=rep(types,boot)) %>%
      group_by(Variant_Effect) %>%
      dplyr::summarise(sd=sd(ratio),
                       mean=mean(ratio),
                       median=median(ratio),
                       upper_conf_95=quantile(ratio, probs = c(0.025,0.975))[2],
                       lower_conf_95=quantile(ratio, probs = c(0.025,0.975))[1],
                       Pop="LHIP")
    bootstrap_lhip$Variant_Effect <- factor(bootstrap_lhip$Variant_Effect,
                                            labels=types,
                                            levels=types)
    
    ratio <- lhip_base %>%
      group_by(Variant_Effect,ROH_50) %>%
      dplyr::summarise(tabulate=n()) %>% 
      spread(key="ROH_50",value="tabulate") %>% 
      mutate(ratio = Lower / (Lower+Upper))
    
    bootstrap_lhip <- merge(ratio,bootstrap_lhip)
    
  }
  
}

bootstraps <- rbind(bootstrap_lhisi,bootstrap_lhip)

bootstraps$Variant_Effect <- factor(bootstraps$Variant_Effect,
                                    levels=c("HIGH","MODERATE","LOW","MODIFIER"),
                                    labels=c("HIGH","MODERATE","LOW","MODIFIER"))


p <- bootstraps %>%
  ggplot() +
  geom_errorbar(aes(y=Pop,x=ratio,xmin=lower_conf_95,xmax=upper_conf_95,fill=Variant_Effect),
                position=position_dodge(width=0.6),width=0,
                filter(bootstraps,Pop=="LHISI")) +
  geom_point(aes(y=Pop,x=ratio,xmin=lower_conf_95,xmax=upper_conf_95,fill=Variant_Effect),
             position=position_dodge(width=0.6),pch=21,size=4.5,
             filter(bootstraps,Pop=="LHISI")) +
  scale_fill_manual(values=rev(new_colours[,3])) +
  new_scale_fill() +
  geom_errorbar(aes(y=Pop,x=ratio,xmin=lower_conf_95,xmax=upper_conf_95,fill=Variant_Effect),
                position=position_dodge(width=0.6),width=0,
                filter(bootstraps,Pop=="LHIP")) +
  geom_point(aes(y=Pop,x=ratio,xmin=lower_conf_95,xmax=upper_conf_95,fill=Variant_Effect),
             position=position_dodge(width=0.6),pch=21,size=4.5,
             filter(bootstraps,Pop=="LHIP")) +
  scale_fill_manual(values=rev(new_colours[,2])) +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.text.x=element_text(size=10),
        axis.text.y=element_blank(),
        legend.position="none") + 
  scale_x_continuous(
    limits=c(0.51,0.78),
    breaks=c(0.55,0.65,0.75),
    labels=c(55,65,75)) +
  labs(x="% variants in\nlow ROH regions") ; p
# Plot

png(paste0(FIGURE_DIR,"/muts_in_roh_",format(Sys.time(),"%Y%m%d"),".png"),
    width=6,height=3,units='in',res=300)
plot(p)
dev.off()
