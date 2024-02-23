### Analysing ngsParalog output

lr <- read.table("lhisi_wgs")
lr$pval <- 0.5*pchisq(lr$LR,df=1,lower.tail=FALSE) # append column of p-values
lr$pval.adj.fdr <- p.adjust(lr$pval, method="fdr") # p-values adjusted for number of tested sites
lr$pval.adj.bon <- p.adjust(lr$pval, method="bonferroni")


# The 7th column of the lr data.frame is the adjusted p-value for rejecting the null hypothesis that reads
# covering the site derive from a single locus. Of course you can use any p-value adjustment of your
# choosing, e.g. "fdr".

# generate list of sites that don't show evidence of mismapping at 0.05 significance level:
qc.sites <- lr[-which(lr$pval.adj < 0.05),1:2]

dim(qc.sites)
dim(lr)


values = lr$pval.adj.fdr
breaks = seq(0, 1, by=0.01) 
values.cut = cut(values, breaks, right=FALSE) 
values.freq = table(values.cut)
values.cumfreq = cumsum(values.freq) / max(cumsum(values.freq))

sum(values > 0.05) / length(values)

library(ggplot2)
library(magrittr)
library(tidyr)
library(dplyr)

HOME_DIR="/Volumes/Alter/LHISI"

### Plotting values along continuous scaffolds

# Read in scaffold information
scaffolds <- read.table(paste0(HOME_DIR,"/References/major_scaffold_regions.bed"))[,c(1,3)]
scaffoldsV3 <- scaffolds$V3 + 1
colnames(scaffolds) <- c("Scaffold","Length")
# Remove sex chromosome
scaffolds <- scaffolds[-4,]

# Get start and end in global coordinates
scaffolds$CumEnd <- cumsum(as.numeric(scaffolds$Length))
scaffolds$CumStart <- c(1,lag(scaffolds$CumEnd)[-1]+1)

# Read in data
lr <- read.table("lhisi_wgs")
colnames(lr) <- c("Scaffold","Position","-log(null)","-log(alternative)","LR")
lr$pval <- 0.5*pchisq(lr$LR,df=1,lower.tail=FALSE) # append column of p-values
lr$pval.adj.fdr <- p.adjust(lr$pval, method="fdr") # p-values adjusted for number of tested sites
lr$pval.adj.bon <- p.adjust(lr$pval, method="bonferroni")

# Give every position a cumulative position
lr <- merge(scaffolds,lr) %>% mutate(CumPosition=Position+CumStart-1)

# Get positions of axis labels
axis_set <- scaffolds %>% 
  group_by(Scaffold) %>% 
  summarize(center = (CumStart+CumEnd)/2)

# Set y-axis minimum and maximum
ymin=0
ymax <- lr %>% 
  mutate(ylim = max(LR) + (max(LR)/100)) %>% 
  summarise(ylim=first(ylim)) %>%
  pull(ylim)

# Data.frame to make panel background rectangles
rect_breaks <- seq(from=ymin,to=ymax,by=(ymax-ymin)/5)
rect <- data.frame(ymin=rect_breaks[c(2,4)],
                   ymax=rect_breaks[c(3,5)],
                   xmax=Inf,xmin=-Inf)

# Make y_labels from this, max and min slightly nudged inward 
y_labels <- rect_breaks
y_labels[1] <- rect_breaks[1] + (rect_breaks[6]/100)
y_labels[6] <- rect_breaks[6] - (rect_breaks[6]/100)

ggplot(lr,aes(x=CumPosition,y=LR,colour=as.factor(Scaffold))) + 
  geom_rect(inherit.aes=NULL,data=rect,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            colour="gray95",fill="gray95") +
  geom_point(alpha=0.75) + 
  scale_x_continuous(label = axis_set$Scaffold, breaks = axis_set$center,
                     limits=c(1-(max(scaffolds$CumEnd)*0.01),
                              max(scaffolds$CumEnd)+(max(scaffolds$CumEnd)*0.01)),
                     expand=c(0,0)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, ymax)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Likelihood ratio\nof mismapping") + 
  theme( 
    legend.position = "none",
    panel.grid=element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )

