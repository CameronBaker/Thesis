#Cameron Baker
#This R script handles the generation of statistical results related
#to the alignment scores for motifs passed in.
#Requires a full run data dump of the motif generator
#Simple pipe is as follows
#GenerateMotif -> EvalMotif -> boxplot_motifs

#Output is a number of boxplots relating to alignment scores
#and the flagging of motif pairs that are not signicantly different into an outfile

boxplot_motifs <- function(filename){
	df <- read.csv(filename,header=FALSE)
	if(length(grep("pf",filename)) > 0){
		desig <- "PF"
	}else{
		desig <- "UP"
	}
	motifs <- unique(df[,1])
	for(motif in motifs){
		single_motif = df[which(df[,1] == motif),c(2,4,5)]
		num_motifs = length(unique(single_motif[,1]))
		title <- gsub("[.]","_",motif)

		other_motifs = unique(single_motif[,1])[which(unique(single_motif[,1]) != motif)]
		for(o_m in other_motifs){
			rmsd_set = single_motif[which(single_motif[,1] == motif),"V4"]
			rmsd_set_other = single_motif[which(single_motif[,1] == o_m),"V4"]
			set_mean = mean(rmsd_set)
			q <- quantile(rmsd_set_other,c(0.25,0.75))
			if(set_mean > q[1] && set_mean < q[2]){
				write(paste(desig,motif,o_m),"tocheck.txt",append=TRUE,sep='\n')
			}
		}	

		jpeg(paste(title,".jpg",sep=""),width=200*num_motifs)
		boxplot(V4~V2,data=single_motif,ylim=c(0,2),main=paste(desig,": ",motif,sep=""),
			  xlab="Protein Designation",ylab="RMSD")
		dev.off()

		jpeg(paste(title,"_super.jpg",sep=""),width=200*num_motifs)
		boxplot(V5~V2,data=single_motif,ylim=c(0,0.5),main=paste(desig," super: ",motif,sep=""),
			  xlab="Protein Designation",ylab="RMSD")
		dev.off()
	}
	print("Done!")
}

file.remove("tocheck.txt")
flist <- list.files()[grep(".txt",list.files())]
for(file in flist){
	boxplot_motifs(file)
}
