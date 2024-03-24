#Rscript  homoout2list.R --freq 2
library("optparse")

option_list = list(
  make_option("--freq", type="character", default=NULL,
              help="kmer to plot", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

myfreq<-as.numeric(as.character(opt$freq))

mytable<-read.table("homo_blastgffgen_out.txt",header=F,sep="\t")
mycount<-as.data.frame(table(mytable$V1))
mykeep<-mycount[which(mycount$Freq>=myfreq),"Var1"]
write.table(mykeep,file="homo_freq_list.txt",quote=F,col.names=F,row.names=F,sep="\t")
