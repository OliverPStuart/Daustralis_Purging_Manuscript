# Source code for colour scheme for WGS paper

colors <- read.delim(paste0(REF_DIR,"/three_pop_colour_map.txt"),
                     header=T,stringsAsFactors=F,sep="\t")
cols <- colors[c(1,2,3),"Set1"]
names(cols) <- c("Wild","LHIP","LHISI")
my_fill_2 <- scale_fill_manual(name = "Population", values = cols[2:3],
                               breaks=c("LHIP","LHISI"),
                               labels=c("Hybrid","Captive"))
my_col_2 <- scale_color_manual(name = "Population", values = cols[2:3],
                               breaks=c("LHIP","LHISI"),
                               labels=c("Hybrid","Captive"))
my_fill_3 <- scale_fill_manual(name = "Population", values = cols,
                               breaks=c("Wild","LHIP","LHISI"),
                               labels=c("Wild","Hybrid","Captive"))
my_col_3 <- scale_color_manual(name = "Population", values = cols,
                               breaks=c("Wild","LHIP","LHISI"),
                               labels=c("Wild","Hybrid","Captive"))
