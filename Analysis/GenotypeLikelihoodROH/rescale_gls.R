
# Library
# Only plyr for mapvalues
suppressWarnings(library(plyr))

# Get input and check

file <- commandArgs(trailingOnly=T)

if(length(file) != 1){
  print("More than one file input, quitting...")
  stop()
} else {
  print(paste0("Converting file ",file))
}

# Get prefix for file writing

prefix=gsub("\\.beagle\\.gz","",file)

# Read table

raw_data <- read.table(file,header=T,check.names=F)

# Get likelihoods

gls <- raw_data[,c(4:ncol(raw_data))]

# First do a matrix calculation to turn all of these
# into Phred-scaled likelihoods

phred_scale <- function(x){
  
  # Add a tiny amount to any 0 values.
  # Cannot take the log of this. Then
  # also round up to integer.
  return(ceiling(-10*log(x+0.000001)))
  
}

pls <- apply(gls,2,phred_scale)

# Then go row by row, column by column and normalise
# so that the most likely has a phred likelihood of 0

n_ind <- dim(pls)[2] / 3
n_sites <- dim(pls)[1]

for(ind in 1:n_ind){
  
  for(site in 1:n_sites){
    
  pls_i <- pls[site,(ind*3-2):(ind*3)]
  pls_i <- pls_i - min(pls_i)
  
  pls[site,(ind*3-2):(ind*3)] <- pls_i

  }
  
}

# This is really fast, so should be tractable for
# larger datasets

# Now reformat the first few columns

head <- data.frame(chr=gsub("_\\d+$","",raw_data$marker),
           id=raw_data$marker,
           pos=gsub("^[^_]*_", "",raw_data$marker),
           allele1=mapvalues(raw_data$allele1,from=c(0,1,2,3),to=c("A","C","G","T")),
           allele2=mapvalues(raw_data$allele2,from=c(0,1,2,3),to=c("A","C","G","T")))

# cbind onto to the PLs

pls <- cbind(head,pls)

# Write the table and gzip it for space

write.table(pls,paste0(prefix,".pl"),quote=F,col.names=F,row.names=F,sep="\t")
system(paste0("gzip ",prefix,".pl"))

