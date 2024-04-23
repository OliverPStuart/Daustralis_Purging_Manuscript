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
library(ggbreak)

options(dplyr.summarise.inform = FALSE)

### Final figures and analyses for deleterious mutations...

### Making data.frame

# Function for complement
`%ni%` <- Negate(`%in%`)

# Read data in, rename for ease of use
Captive <- read.table("../AlleleFrequencies/lhisi.mafs.gz",header=T,stringsAsFactors=F,sep="\t")[,c(1,2,6,7)]
Hybrid <- read.table("../AlleleFrequencies/lhip.mafs.gz",header=T,stringsAsFactors=F,sep="\t")[,c(1,2,6,7)]
Wild <- read.table("../AlleleFrequencies/Wild.mafs.gz",header=T,stringsAsFactors=F,sep="\t")[,c(1,2,6,7)]
colnames(Captive) <- c("Scaffold","Position","MAF_Captive","N_Captive")
colnames(Hybrid) <- c("Scaffold","Position","MAF_Hybrid","N_Hybrid")
colnames(Wild) <- c("Scaffold","Position","MAF_Wild","N_Wild")

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

Wild$fate_Wild <- sapply(Wild$MAF_Wild,classify_maf)
Captive$fate_Captive <- sapply(Captive$MAF_Captive,classify_maf)
Hybrid$fate_Hybrid <- sapply(Hybrid$MAF_Hybrid,classify_maf)

# Now get all of these into the effects dataset
effects <- merge(effects,Wild,all.x=T)
effects <- merge(effects,Captive,all.x=T)
effects <- merge(effects,Hybrid,all.x=T)

# Make an ID column for later
effects$ID <- paste0(effects$Scaffold,"_",effects$Position)

rm(Captive,Hybrid,Wild)

effects <- effects %>% mutate(captive_fate=ifelse((fate_Wild == "absent" & 
                                                     fate_Captive != "absent" &  
                                                     fate_Hybrid != "absent"), "common",
                                                  ifelse((fate_Wild == "absent" &  
                                                            fate_Captive != "absent" &  
                                                            fate_Hybrid == "absent"), "Captive",
                                                         ifelse((fate_Wild == "absent" &  
                                                                   fate_Captive == "absent" &  
                                                                   fate_Hybrid != "absent"),  "Hybrid",NA))))

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


# Classifying variants based on presence in putative ROH in the three populations

Hybrid <- read.table("LHIP.cov",header=F,stringsAsFactors=F,sep="\t")
Captive <- read.table("LHISI.cov",header=F,stringsAsFactors=F,sep="\t")
Wild <- read.table("Wild.cov",header=F,stringsAsFactors=F,sep="\t")
colnames(Hybrid) <- c("Scaffold","Position","Coverage_Hybrid")
colnames(Captive) <- c("Scaffold","Position","Coverage_Captive")
colnames(Wild) <- c("Scaffold","Position","Coverage_Wild")

tmp <- merge(Hybrid,Captive)
tmp <- merge(tmp,Wild)

effects <- merge(effects,tmp,all.x=T)

rm(tmp,Captive,Hybrid,Wild)

# Some refactorisation to order specific variables

effects$Variant_Effect <- factor(effects$Variant_Effect,
                                 levels=c("MODIFIER","LOW","MODERATE","HIGH"),
                                 labels=c("MODIFIER","LOW","MODERATE","HIGH"))

# Colour map

# Get the file of colour schemes
colors <- read.delim(paste0(REF_DIR,"/three_pop_colour_map.txt"),
                     header=T,stringsAsFactors=F,sep="\t")
cols <- colors[c(1,2,3),"Set1"]
names(cols) <- c("Wild","Hybrid","Captive")
my_fill_2 <- scale_fill_manual(name = "Population", values = cols[2:3])
my_col_2 <- scale_color_manual(name = "Population", values = cols[2:3])
my_fill_3 <- scale_fill_manual(name = "Population", values = cols)
my_col_3 <- scale_color_manual(name = "Population", values = cols)

## Change AED threshold
#
#effects <- effects[effects$AED < 0.7 | is.na(effects$AED),]

# Make figure directory
if(!file.exists("figures")){
  dir.create("figures")
}

setwd("figures")

# Plot of proportion of different effect classes segregating in three populations
wild <- effects %>% filter(fate_Wild == "segregating") %>% pull(Variant_Effect)
hybrid <- effects %>% filter(fate_Hybrid == "segregating") %>% pull(Variant_Effect)
captive <- effects %>% filter(fate_Captive == "segregating") %>% pull(Variant_Effect)

all <- data.frame(effect = c(wild,hybrid,captive),
                  pop = c(rep("Wild",length(wild)),
                          rep("Hybrid",length(hybrid)),
                          rep("Captive",length(captive))))
all$pop <- factor(all$pop,levels=c("Wild","Hybrid","Captive"),
                  labels=c("Wild","Hybrid","Captive"))

png("variant_effect_breakdown_by_pop.png",res=330,width=8,height=4,units='in')
all %>% ggplot(aes(y=pop,fill=effect)) + 
  geom_bar(position="fill",colour="black",width=0.95) + 
  scale_x_break(c(0.025,0.95),scales=0.2,ticklabels=c(0.95,1.00),expand=c(0,0),space=0.7) + 
  scale_y_discrete(expand=c(0,0)) +
  theme_bw() + 
  theme(panel.border = element_rect(colour="black"),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        panel.grid=element_blank(),
        axis.text.x=element_text(angle=30,hjust=1,vjust=1),
        axis.text.y=element_text(size=12),
        axis.text.y.right = element_blank()) + 
  scale_fill_brewer(palette = "YlGn",name="Variant\neffect",
                    breaks=c("MODIFIER","LOW","MODERATE","HIGH"),
                    labels=c("NON-CODING","LOW","MODERATE","HIGH"),
  )
dev.off()

# Frequency difference of different mutation classes in captivity compared to the wild

Wild <- effects %>% filter(fate_Wild != "absent")

# Calculate difference in frequency
Wild$Captive_change <- Wild$MAF_Captive - Wild$MAF_Wild
Wild$Hybrid_change <- Wild$MAF_Hybrid - Wild$MAF_Wild
Wild <- Wild %>% drop_na(Captive_change,Hybrid_change)

# Plotting mean+sd difference per variant effect class
# Excluding those variants in each population which are absent
tmp_capt <- Wild %>% group_by(Variant_Effect) %>% 
  filter(fate_Captive != "absent") %>%
  dplyr::summarise(mean_Captive=mean(Captive_change),
                   sd_Captive=sd(Captive_change)) %>%
  gather(key, value, -Variant_Effect) %>%
  extract(key, c("stat", "pop"), "(.+)_(.+)") %>%
  spread(stat, value)
tmp_hybr <- Wild %>% group_by(Variant_Effect) %>% 
  filter(fate_Hybrid != "absent") %>%
  dplyr::summarise(mean_Hybrid=mean(Hybrid_change),
                   sd_Hybrid=sd(Hybrid_change)) %>%
  gather(key, value, -Variant_Effect) %>%
  extract(key, c("stat", "pop"), "(.+)_(.+)") %>%
  spread(stat, value)

change_data <- rbind(tmp_capt,tmp_hybr)

p <- change_data %>% ggplot(aes(x=Variant_Effect,y=mean,fill=pop)) + 
  scale_y_continuous(limits=c(-0.7,0.1)) + 
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0,size=0.7,position=position_dodge(width=0.3)) + 
  geom_line(aes(group=pop),linetype="dashed") +
  geom_point(size=4.5,pch=21,position=position_dodge(width=0.3)) + 
  geom_hline(yintercept=0,linetype="dashed",colour="black",size=1) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  labs(x="Variant\neffect",y="Frequency difference\nin captivity") + 
  my_fill_2 + 
  scale_x_discrete(expand=c(-0.3,1.4))

png("frequency_difference_by_pop.png",res=300,width=6,height=7,units='in')
plot(p)
dev.off()

p <- change_data %>% ggplot(aes(fill=Variant_Effect,y=mean,x=pop)) + 
  scale_y_continuous(limits=c(-0.7,0.1)) + 
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0,size=0.7,position=position_dodge(width=0.6)) + 
  geom_point(size=4.5,pch=21,position=position_dodge(width=0.6)) + 
  geom_hline(yintercept=0,linetype="dashed",colour="black",size=1) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  labs(x="Population",y="Frequency difference\nin captivity") + 
  scale_fill_brewer(palette = "YlGn",name="Variant\neffect")

png("frequency_difference_by_type.png",res=300,width=6,height=7,units='in')
plot(p)
dev.off()

p <- change_data %>% ggplot(aes(fill=factor(Variant_Effect,levels=c("HIGH","MODERATE","LOW","MODIFIER")),
                           x=mean,y=pop)) + 
  scale_x_continuous(limits=c(-0.7,0.1)) + 
  geom_errorbar(aes(xmin=mean-sd,xmax=mean+sd),width=0,size=0.7,position=position_dodge(width=0.6)) + 
  geom_point(size=4.5,pch=21,position=position_dodge(width=0.6)) + 
  geom_vline(xintercept=0,linetype="dashed",colour="black",size=1) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.minor=element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(size=10),
        panel.grid=element_blank()) + 
  labs(x="Frequency difference\nin captivity") + 
  scale_fill_brewer(palette = "YlGn",name="Variant\neffect",
                    breaks=c("MODIFIER","LOW","MODERATE","HIGH"),
                    labels=c("NON-CODING","LOW","MODERATE","HIGH"),
#                    guide = guide_legend(reverse = TRUE)
                    direction=-1)

png("frequency_difference_by_type_flipped.png",res=300,width=7,height=3,units='in')
plot(p)
dev.off()

# Ratio of deleterious variant presence inside/outside of high ROH regions
# Need to standardise this by presence of ROH in the first place....
# Would that just alter the difference between populations and not the populations themselves?

capt <- effects %>% filter(
  fate_Captive == "segregating",
  #Variant_Effect != "MODIFIER"
) %>%
  mutate(ROH_50 = ifelse(Coverage_Captive > 0.5,"Upper","Lower")) %>%
  group_by(Variant_Effect,ROH_50) %>%
  dplyr::summarise(tabulate=n()) %>% 
  spread(key="ROH_50",value="tabulate") %>% 
  mutate(ratio = Lower / Upper,pop="Captive")

hybr <- effects %>% filter(
  fate_Hybrid == "segregating",
  #Variant_Effect != "MODIFIER"
) %>%
  mutate(ROH_50 = ifelse(Coverage_Hybrid > 0.5,"Upper","Lower")) %>%
  group_by(Variant_Effect,ROH_50) %>%
  dplyr::summarise(tabulate=n()) %>% 
  spread(key="ROH_50",value="tabulate") %>% 
  mutate(ratio = Lower / Upper,pop="Hybrid")

p <- rbind(capt,hybr) %>%
  ggplot(aes(y=pop,x=ratio,fill=factor(Variant_Effect,levels=c("HIGH","MODERATE","LOW","MODIFIER")))) + 
  geom_linerange( aes(y=pop, xmin=1, xmax=ratio),position=position_dodge(width=0.5)) +
  geom_point(size=4.5,pch=21,position=position_dodge(width=0.5)) + 
  scale_fill_brewer(palette = "YlGn",name="Variant\neffect",
                    breaks=c("MODIFIER","LOW","MODERATE","HIGH"),
                    labels=c("NON-CODING","LOW","MODERATE","HIGH"),
                    direction=-1) + 
  scale_x_continuous(limits=c(1,2.5),expand=c(0,0)) + 
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=10)) + 
  labs(x="Depletion of variants\ninside ROH")

png("muts_in_roh_ratio.png",res=300,width=7,height=3,units='in')
plot(p)
dev.off()


# Getting some sampling error intervals for these with bootstrapping

capt_base <- effects %>% filter(
  fate_Captive == "segregating",
  #Variant_Effect != "MODIFIER"
) %>%
  mutate(ROH_50 = ifelse(Coverage_Captive > 0.5,"Upper","Lower"))

hybr_base <- effects %>% filter(
  fate_Hybrid == "segregating"
) %>%
  mutate(ROH_50 = ifelse(Coverage_Hybrid > 0.5,"Upper","Lower"))

# Now bootstrap based on proportions of variant classes

# Number of bootstraps
boot=100
# Get variant vector
types <- c("MODIFIER","LOW","MODERATE","HIGH")

# Captive
vec <- c()
for(i in 1:boot){

  for(j in 1:length(types)){
    
    tmp <- capt_base[capt_base$Variant_Effect == types[j],]
      
    val <- tmp[sample(x=1:nrow(tmp),size=nrow(tmp),replace=T),] %>%
      group_by(ROH_50) %>%
      dplyr::summarise(tabulate=n()) %>% 
      spread(key="ROH_50",value="tabulate") %>% 
      mutate(ratio = Lower / (Lower+Upper)) %>%
      pull(ratio)
  
    vec <- c(vec,val)
  
  }
  
  if(i == boot){
    bootstrap_capt <- data.frame(ratio=vec,
                                 Variant_Effect=rep(types,boot)) %>%
      group_by(Variant_Effect) %>%
      dplyr::summarise(sd=sd(ratio),
                       mean=mean(ratio),
                       median=median(ratio),
                       upper_conf_95=quantile(ratio, probs = c(0.025,0.975))[2],
                       lower_conf_95=quantile(ratio, probs = c(0.025,0.975))[1],
                       Pop="Captive")
    
    bootstrap_capt$Variant_Effect <- factor(bootstrap_capt$Variant_Effect,
                                            labels=types,
                                            levels=types)
    
    ratio <- capt_base %>%
      group_by(Variant_Effect,ROH_50) %>%
      dplyr::summarise(tabulate=n()) %>% 
      spread(key="ROH_50",value="tabulate") %>% 
      mutate(ratio = Lower / (Lower+Upper))
    
    bootstrap_capt <- merge(ratio,bootstrap_capt)
    
  }
  
}

# Hybrid
vec <- c()
for(i in 1:boot){
  
  for(j in 1:length(types)){
    
    tmp <- hybr_base[hybr_base$Variant_Effect == types[j],]
    
    val <- tmp[sample(x=1:nrow(tmp),size=nrow(tmp),replace=T),] %>%
      group_by(ROH_50) %>%
      dplyr::summarise(tabulate=n()) %>% 
      spread(key="ROH_50",value="tabulate") %>% 
      mutate(ratio = Lower / (Lower+Upper)) %>%
      pull(ratio)
    
    vec <- c(vec,val)
    
  }
  
  if(i == boot){
    bootstrap_hybr <- data.frame(ratio=vec,
                                 Variant_Effect=rep(types,boot)) %>%
      group_by(Variant_Effect) %>%
      dplyr::summarise(sd=sd(ratio),
                       mean=mean(ratio),
                       median=median(ratio),
                       upper_conf_95=quantile(ratio, probs = c(0.025,0.975))[2],
                       lower_conf_95=quantile(ratio, probs = c(0.025,0.975))[1],
                       Pop="Hybrid")
    bootstrap_hybr$Variant_Effect <- factor(bootstrap_hybr$Variant_Effect,
                                            labels=types,
                                            levels=types)
    
    ratio <- hybr_base %>%
      group_by(Variant_Effect,ROH_50) %>%
      dplyr::summarise(tabulate=n()) %>% 
      spread(key="ROH_50",value="tabulate") %>% 
      mutate(ratio = Lower / (Lower+Upper))
    
    bootstrap_hybr <- merge(ratio,bootstrap_hybr)
    
  }
  
}

##### WORKING ON COLOUR SCHEME

rbind(bootstrap_capt,bootstrap_hybr) %>%
  ggplot(aes(y=Pop,x=ratio,xmin=lower_conf_95,xmax=upper_conf_95,shape=Variant_Effect,fill=Pop)) +
  geom_errorbar(position=position_dodge(width=0.5),width=0) +
  geom_point(position=position_dodge(width=0.5),size=3) +
  my_fill_2 +
  scale_shape_manual(values=21:24) +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=10)) + 
  labs(x="% variants in\nlow ROH regions") + 
  scale_x_continuous(limits=c(0.5,0.8),expand=c(0,0))

# Two dimensional colour scheme?

# Now what if we combine all coding variant types
# This should improve statistical power while eroding some subtlety

capt_base <- effects %>% filter(
  fate_Captive == "segregating",
  #Variant_Effect != "MODIFIER"
) %>%
  mutate(ROH_50 = ifelse(Coverage_Captive > 0.5,"Upper","Lower"),
         coding = ifelse(Variant_Effect == "MODIFIER","NEUTRAL","FUNCTIONAL"))

hybr_base <- effects %>% filter(
  fate_Hybrid == "segregating"
) %>%
  mutate(ROH_50 = ifelse(Coverage_Hybrid > 0.5,"Upper","Lower"),
         coding = ifelse(Variant_Effect == "MODIFIER","NEUTRAL","FUNCTIONAL"))

# Number of bootstraps
boot=100
# Get variant vector
types <- c("NEUTRAL","FUNCTIONAL")

# Captive
vec <- c()
for(i in 1:boot){
  
  for(j in 1:length(types)){
    
    tmp <- capt_base[capt_base$coding == types[j],]
    
    val <- tmp[sample(x=1:nrow(tmp),size=nrow(tmp),replace=T),] %>%
      group_by(ROH_50) %>%
      dplyr::summarise(tabulate=n()) %>% 
      spread(key="ROH_50",value="tabulate") %>% 
      mutate(ratio = Lower / (Lower+Upper)) %>%
      pull(ratio)
    
    vec <- c(vec,val)
    
  }
  
  if(i == boot){
    bootstrap_capt <- data.frame(ratio=vec,
                                 coding=rep(types,boot)) %>%
      group_by(coding) %>%
      dplyr::summarise(sd=sd(ratio),
                       mean=mean(ratio),
                       median=median(ratio),
                       upper_conf_95=quantile(ratio, probs = c(0.025,0.975))[2],
                       lower_conf_95=quantile(ratio, probs = c(0.025,0.975))[1],
                       Pop="Captive",
                       min=min(ratio),
                       max=max(ratio))
    
    bootstrap_capt$coding <- factor(bootstrap_capt$coding,
                                            labels=types,
                                            levels=types)
    
    ratio <- capt_base %>%
      group_by(coding,ROH_50) %>%
      dplyr::summarise(tabulate=n()) %>% 
      spread(key="ROH_50",value="tabulate") %>% 
      mutate(ratio = Lower / (Lower+Upper))
    
    bootstrap_capt <- merge(ratio,bootstrap_capt)
    
  }
  
}

# Hybrid
vec <- c()
for(i in 1:boot){
  
  for(j in 1:length(types)){
    
    tmp <- hybr_base[hybr_base$coding == types[j],]
    
    val <- tmp[sample(x=1:nrow(tmp),size=nrow(tmp),replace=T),] %>%
      group_by(ROH_50) %>%
      dplyr::summarise(tabulate=n()) %>% 
      spread(key="ROH_50",value="tabulate") %>% 
      mutate(ratio = Lower / (Lower+Upper)) %>%
      pull(ratio)
    
    vec <- c(vec,val)
    
  }
  
  if(i == boot){
    bootstrap_hybr <- data.frame(ratio=vec,
                                 coding=rep(types,boot)) %>%
      group_by(coding) %>%
      dplyr::summarise(sd=sd(ratio),
                       mean=mean(ratio),
                       median=median(ratio),
                       upper_conf_95=quantile(ratio, probs = c(0.025,0.975))[2],
                       lower_conf_95=quantile(ratio, probs = c(0.025,0.975))[1],
                       Pop="Hybrid",
                       min=min(ratio),
                       max=max(ratio))
    bootstrap_hybr$coding <- factor(bootstrap_hybr$coding,
                                            labels=types,
                                            levels=types)
    
    ratio <- hybr_base %>%
      group_by(coding,ROH_50) %>%
      dplyr::summarise(tabulate=n()) %>% 
      spread(key="ROH_50",value="tabulate") %>% 
      mutate(ratio = Lower / (Lower+Upper))
    
    bootstrap_hybr <- merge(ratio,bootstrap_hybr)
    
  }
  
}

rbind(bootstrap_capt,bootstrap_hybr) %>%
  ggplot(aes(y=Pop,x=ratio,xmin=lower_conf_95,xmax=upper_conf_95,fill=coding)) +
  geom_errorbar(position=position_dodge(width=0.5),width=0) +
  geom_point(position=position_dodge(width=0.5),pch=21,size=3) +
  scale_fill_brewer(palette = "YlGn") +
  theme_bw() + 
  theme(axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        axis.text=element_text(size=10)) + 
  labs(x="% variants in\nlow ROH regions")


