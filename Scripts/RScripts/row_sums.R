d <- read.table("tmp")
means <- rowSums(d) / ncol(d)
writeLines(con="means",text=as.character(means))
