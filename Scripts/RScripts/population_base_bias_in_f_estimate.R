### Environment

# Directory

HOME_DIR="/Volumes/Alter/LHISI"
REF_DIR=paste0(HOME_DIR,"/References")

#HOME_DIR="/home/ollie/LHISI_Thesis_Work"
#REF_DIR=paste0("~/References")

WORKING_DIR=paste0(HOME_DIR,"/Analyses/GenotypeLikelihoodROH")

setwd(WORKING_DIR)

# Libraries
library(RZooRoH)
library(tidyr)
library(ggplot2)
library(plyr)
library(dplyr)

# Read in all data objects

all <- load("all_results/all_model_results.RData")
lhip <- load("lhip_results/lhip_model_results.RData")
lhisi <- load("lhisi_results/lhisi_model_results.RData")

# Compare parameters per population on F

# Empty data.frame

output <- data.frame(fs=numeric(),names=character(),
           k=numeric(),rate=numeric(),pop=character())

pops <- c("all","lhisi","lhip")
rates <- c(2,5,10)
classes <- c(6:10)

for(i in 1:length(pops)){
  for(j in 1:length(rates)){
    for(class in 1:length(classes)){
  
      
      pop <- pops[i]
      rate <- rates[j]
      k <- classes[class]
    
  fs <- eval(parse(text=paste0("1- ",pop,"_K",k,"_rate",rate,"_results@realized[,",k,"]")))
  names <- readLines(paste0(pop,"_bam_list.txt"))
  names <- data.frame(Sample=gsub("\\.bam","",gsub(".*/","",names)),
                      id=1:length(names))
  
  output <- rbind(output,
        data.frame(fs=fs,names=names$Sample,
             k=k,rate=rate,pop=pop))
  
  }
  }
}
  

names_lhip <- readLines(paste0("lhip_bam_list.txt"))
names_lhip <- gsub("\\.bam","",gsub(".*/","",names_lhip))
names_lhisi <- readLines(paste0("lhisi_bam_list.txt"))
names_lhisi <- gsub("\\.bam","",gsub(".*/","",names_lhisi))

# Now plot F for lhisi,lhip individuals
# Comparing all or population specific allele frequencies

output %>% dplyr::filter(names %in% names_lhip) %>%
  dplyr::filter(rate==2) %>%
  ggplot(aes(x=k,col=names,y=fs,group=names,shape=pop)) + 
  geom_point(position=position_dodge(width=0.2)) + 
  scale_y_continuous(limits=c(0,1))

output %>% dplyr::filter(names %in% names_lhisi) %>%
  dplyr::filter(rate==2) %>%
  ggplot(aes(x=k,col=names,y=fs,group=names,shape=pop)) + 
  geom_point(position=position_dodge(width=0.2)) + 
  scale_y_continuous(limits=c(0,1))

# K has no effect on the values
# Rates has no effect on the values
# For all future comparisons, we are just looking at K=9, rate=2

# Calculate a mean difference in f per individual

if(!file.exists("figures")){
  dir.create("figures")
  }

setwd("figures")
  
p <- output %>% dplyr::filter(names %in% names_lhip) %>%
  pivot_wider(names_from=pop,values_from=fs) %>%
  mutate(difference = all-lhip,rate_c=as.character(rate),k_c=as.character(k)) %>%
  ggplot(aes(x=rate_c,y=difference,fill=k_c)) + 
  geom_boxplot() + 
  scale_y_continuous(limits=c(-0.15,0.15)) + 
  geom_hline(yintercept=0,linetype="dashed") +
  labs(title = "Excess F in LHIP when\nusing all as reference",
       x = "Stopping rate",
       y = "Difference",
       fill = "K")

png("lhip_f_bias.png",res=300,width=6,height=6,units='in')
plot(p)
dev.off()

p <- output %>% dplyr::filter(names %in% names_lhisi) %>%
  pivot_wider(names_from=pop,values_from=fs) %>%
  mutate(difference = all-lhisi,rate_c=as.character(rate),k_c=as.character(k)) %>%
  ggplot(aes(x=rate_c,y=difference,fill=k_c)) + 
  geom_boxplot() + 
  scale_y_continuous(limits=c(-0.15,0.15)) + 
  geom_hline(yintercept=0,linetype="dashed")  +
  labs(title = "Excess F in LHISI when\nusing all as reference",
       x = "Stopping rate",
       y = "Difference",
       fill = "K")

png("lhisi_f_bias.png",res=300,width=6,height=6,units='in')
plot(p)
dev.off()

# What is the difference in F between lhip
# and lhisi when only considering population
# specific allele frequencies

details <- read.table(paste0(REF_DIR,"/wgs_sample_details.txt"),
                      header=T,stringsAsFactors = F)
colnames(details) <- c("names","sex","population")
output <- merge(output,details)

output %>% dplyr::filter(names != "VAN2", names != "PAUL4",
                  rate==2,k==9) %>%
  ggplot(aes(x=population,y=fs)) + 
  geom_boxplot() +
  facet_grid(.~pop,space="free",scales="free_x")

# Does this change the statistical significance
tmp1 <- output %>% dplyr::filter(names != "VAN2", names != "PAUL4",
                  rate==2,k==9,population=="LHIP",pop=="lhip") %>%
  pull(fs)
tmp2 <- output %>% dplyr::filter(names != "VAN2", names != "PAUL4",
                          rate==2,k==9,population=="LHISI",pop=="lhisi") %>%
  pull(fs)
t.test(x=tmp1,y=tmp2)


tmp1 <- output %>% dplyr::filter(names != "VAN2", names != "PAUL4", names != "C01220",
                          rate==2,k==9,population=="LHIP",pop=="lhip") %>%
  pull(fs)
tmp2 <- output %>% dplyr::filter(names != "VAN2", names != "PAUL4", names != "C01220",
                          rate==2,k==9,population=="LHISI",pop=="lhisi") %>%
  pull(fs)

t.test(x=tmp1,y=tmp2)

# Nope outcome is the same
# When excluding C01220, difference is significant
