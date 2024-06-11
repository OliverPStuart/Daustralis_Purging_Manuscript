### Environment

# Directory

HOME_DIR="/mnt/data/dayhoff/home/scratch/groups/mikheyev/LHISI"
WORKING_DIR=paste0(HOME_DIR,"/Analysis/GenotypeLikelihoodROH")
REF_DIR=paste0(HOME_DIR,"/References")

setwd(WORKING_DIR)

# Libraries
library(RZooRoH)
library(tidyr)
library(ggplot2)
library(plyr)

# File name

file = commandArgs(trailingOnly=TRUE)
#file="all_scaffold17.pl.gz"

# Get prefix for file writing
# Also get pop (and scaffold?)

prefix=gsub("\\.beagle\\.gz","",file)
pop=strsplit(prefix,split = "_")[[1]][1]

# Get names

names <- readLines(paste0("../../",pop,"_bam_list.txt"))
names <- data.frame(Sample=gsub("\\.bam","",gsub(".*/","",names)),
                    id=1:length(names))
details <- merge(names,read.table(paste0(REF_DIR,"/wgs_sample_details.txt"),
                                  header=T,stringsAsFactors=F,sep="\t"))

# Get plotting order for HBD segments
details$Pop <- mapvalues(details$Pop,from=c("LHIP","LHISI","WILD"),
                         to=c("Hybrid","Captive","Wild"))
details$Pop <- factor(details$Pop,levels=c("Wild","Hybrid","Captive"))
details <- details[order(details$id),]

# Read data

data_gl <- zoodata(genofile = file, zformat = "gl")

# Make empty data frame to attach model comparison
# results to.

models <- data.frame(Sample=character(),
                     Pop=character(),
                     id=character(),
                     K=numeric(),
                     base_rate=numeric(),
                     BIC=numeric(),
                     L=numeric(),
                     F_ROH=numeric(),
                     Pop=character())

# Define vectors of classes and base_rates to try

Ks = c(6,7,8,9,10)
base_rates = c(2,5,10)

for( i in 1:length(Ks) ){
  for ( j in 1:length(base_rates) ){
    
    # Define model
    
    mix1OR <- zoomodel(predefined = T,
                       K = Ks[i],
                       base_rate = base_rates[j])
    
    # Run model
    # Uncomment nT argument if running on cluster

    results <- zoorun(mix1OR,data_gl,
                      method="estem",convem=1e-12,
                      fb = T,minr=1,
                      maxiter=5000,
                      nT=48
                      )
    
    # Attach model results to the data.frame
    
    models <- rbind(models,
                    data.frame(Sample=details$Sample,
                               id=details$id,
                               Pop=details$Pop,
                               K=Ks[i],
                               base_rate=base_rates[j],
                               BIC=results@modbic,
                               L=results@modlik,
                               F_ROH=1-results@realized[,Ks[i]]))
    
    # Store results in new object
    
    assign(paste0(pop,"_K",Ks[i],"_rate",base_rates[j],"_results"),results)
    
  }
}

# Make plots of model results comparing K and rates

dir.create(paste0(pop,"_results"))
setwd(paste0(pop,"_results"))

for(i in 1:length(base_rates)){

p <- ggplot(models[models$base_rate==base_rates[i],],aes(x=K,y=BIC,group=base_rate,color=Pop)) + 
  geom_line() + facet_wrap(~Sample,scales="free_y") + geom_point()
png(paste0(pop,"_rate",base_rates[i],"_model_comparison_BIC.png"),res=300,width=10,height=10,units='in')
plot(p)
dev.off()

p <- ggplot(models[models$base_rate==base_rates[i],],aes(x=K,y=L,group=base_rate,color=Pop)) + 
  geom_line() + facet_wrap(~Sample,scales="free_y") + geom_point()
png(paste0(pop,"_rate",base_rates[i],"_model_comparison_L.png"),res=300,width=10,height=10,units='in')
plot(p)
dev.off()

p <- ggplot(models[models$base_rate==base_rates[i],],aes(x=K,y=F_ROH,group=base_rate,color=Pop)) + 
  geom_line() + facet_wrap(~Sample) + geom_point() +
  scale_y_continuous(limits=c(0,1))
png(paste0(pop,"_rate",base_rates[i],"_model_comparison_F.png"),res=300,width=10,height=10,units='in')
plot(p)
dev.off()

}

# Save .Rdata object for later access to all model results

save(file = paste0(pop,"_model_results.RData"),
     list = ls()[grep(".*rate.*results",ls())])

# Move up again....... just in case

setwd("..")
