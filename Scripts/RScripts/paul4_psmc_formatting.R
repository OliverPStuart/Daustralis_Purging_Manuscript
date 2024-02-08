
args <- commandArgs(trailingOnly = TRUE)
file = args[1]
mu = args[2]
s = args[3]
g = args[4]
i = args[5]

psmc.result<-function(file,i.iteration=25,mu=mu,s=s,g=g)
{
	X<-scan(file=file,what="",sep="\n",quiet=TRUE)
	
	START<-grep("^RD",X)
	END<-grep("^//",X)
	
	X<-X[START[i.iteration+1]:END[i.iteration+1]]
	
	TR<-grep("^TR",X,value=TRUE)
	RS<-grep("^RS",X,value=TRUE)
	
	write(TR,paste("temp.psmc.result_",i))
	theta0<-as.numeric(read.table(paste("temp.psmc.result_",i))[1,2])
	N0<-theta0/4/as.numeric(mu)/as.numeric(s)
	
	write(RS,paste("temp.psmc.result_",i))
	a<-read.table(paste("temp.psmc.result_",i))
	Generation<-as.numeric(2*N0*a[,3])
	Ne<-as.numeric(N0*a[,4])
	
	file.remove(paste("temp.psmc.result_",i))
	
	n.points<-length(Ne)
	YearsAgo<-c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
		Generation[n.points])*as.numeric(g)
	Ne<-c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
		Ne[n.points])
	
	return(data.frame(YearsAgo,Ne,Rep=i))
}

write.table(psmc.result(file=file,mu=mu,s=s,g=g),
            paste0("psmc_formatted_",i,".txt"),
                   quote=F,sep="\t",row.names=F,col.names=F)