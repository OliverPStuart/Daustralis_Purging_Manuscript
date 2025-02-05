### Analysing ngsParalog output

library(ggplot2)
library(magrittr)
library(tidyr)
library(dplyr)

setwd("/Volumes/Alter/Daus_WGS_Paper/Analysis/ngsParalogAnalysis")

### Plotting values along continuous scaffolds

# Read in scaffold information
scaffolds <- read.table("../../References/autosome_regions.bed")[,c(1,3)]
scaffolds$V3 <- scaffolds$V3 + 1
colnames(scaffolds) <- c("Scaffold","Length")

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
  geom_point(alpha=0.5) + 
  scale_x_continuous(label = axis_set$Scaffold, breaks = axis_set$center,
                     limits=c(1-(max(scaffolds$CumEnd)*0.01),
                              max(scaffolds$CumEnd)+(max(scaffolds$CumEnd)*0.01)),
                     expand=c(0,0)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, ymax)) +
  scale_color_manual(values = rep(c("#276FBF", "#188059"), unique(length(axis_set$Scaffold)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Likelihood ratio\nof mismapping") + 
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.grid=element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )

