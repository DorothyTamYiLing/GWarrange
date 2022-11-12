#Usage: Rscript addIS_merge.R --input /home/ubuntu/Dorothy/NCBI_USA_BP_genomes/ISreplacement/IS_replaced_649_genomes/667_IScoor_merge5_ext5.txt --freq 2000 --merge 2000 --ISsize 1000 --gen.size /home/ubuntu/Dorothy/NCBI_USA_BP_genomes/original_USAgenomes/667_genomesize.txt
#merge must be smaller than freq
#size must be smaller than freq
#merge must be smaller than (freq-size)

library("optparse")

option_list = list(
  make_option("--input", type="character", default=NULL, 
              help="input is output of merge_IS_5Oct22.R", metavar="character"),
  make_option("--freq", type="character", default=NULL, 
              help="creating artifical IS start-end coordinates (1000bp) every x bp along the genome", metavar="character"),
  make_option("--merge", type="character", default=NULL, 
              help="maximum number of bp between two adjacent IS for them to be merged into one, this number can be the length of kmer used in GWAS", metavar="character"),
  make_option("--ISsize", type="character", default=NULL, 
              help="size of replacement", metavar="character"),
  make_option("--gen.size", type="character", default=NULL, 
              help="size of genomes, file with no header", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); 

myfreq<-as.numeric(opt$freq)
#myfreq=5000
mymerge<-as.numeric(opt$merge)
#mymerge<-3000
myISsize<-as.numeric(opt$ISsize)
#myISsize<-3000

mygensize<-read.table(opt$gen.size,header=F,sep="\t")
#mygensize<-read.table("/home/ubuntu/Dorothy/NCBI_USA_BP_genomes/original_USAgenomes/667_genomesize.txt",header=F,sep="\t")
myinput<-read.delim(opt$input,header=T,sep="\t")
#myinput<-read.delim("/home/ubuntu/Dorothy/NCBI_USA_BP_genomes/ISreplacement/IS_replaced_649_genomes/667_IScoor_merge5_ext5.txt",header=T,sep="\t")

#creating artifical IS start-end coordinates every myfreq bp along the genome
for(i in 1:nrow(mygensize)){
  mysample<-as.character(mygensize[i,1])
  print(mysample)
  mygenlen<-mygensize[i,2]
  #create the IS midpoint position
  mymidpt<-seq(((myISsize/2)+1000),(mygenlen-(myISsize/2)-1000),by=myfreq)
  #create the start and end matrix for the sample
  mymakeIS<-matrix(0,length(mymidpt),3)
  colnames(mymakeIS)<-c("sseqid","mystart","myend")
  mymakeIS[,1]<-mysample
  mymakeIS[,2]<-as.numeric(mymidpt-(myISsize/2))
  mymakeIS[,3]<-as.numeric(mymidpt+(myISsize/2))
  myinput<-rbind(myinput,mymakeIS)  #add to the original IS coordinates
}

#make the coordinate column numeric
myinput$mystart<-as.numeric(myinput$mystart)
myinput$myend<-as.numeric(myinput$myend)

#now merge again
mymerge_all<-matrix(0,1,3)
colnames(mymerge_all)<-c("sseqid","mystart","myend")

#loop by samples
mylist<-unique(myinput$sseqid)

for(k in 1:length(mylist)){
  print(mylist[k])
  mytable<-myinput[which(myinput$sseqid==mylist[k]),]
  
  #mytable<-myinput[which(myinput$sseqid=="J066"),]
  
  #before running the loop, make sure mystart always < myend, and mystart column is ordered
  mytable<-mytable[order(mytable$mystart),]
  
  #make the merge output matrix   (when the IS is <=10bp apart they are combined into one)
  mymerge_row<-matrix(0,1,2)
  
  mystart_final<-mytable[1,"mystart"]
  myend_final<-mytable[1,"myend"]
  
  for(i in 2:nrow(mytable)){
    mystart<-mytable[i,"mystart"]
    myend<-mytable[i,"myend"]
    if(mystart>=mystart_final & myend<=(myend_final)){ #there should be no "+x" here, but keep it for replicating result
      #print("new row within the range, do not do anything")
    }
    if(mystart<=(myend_final+mymerge) & myend>(myend_final)){ #there should be no "+x" here, but keep it for replicating result
      myend_final<-myend
      #print("end exceed the range but start still within the range, update end only without starting new range")
    }
    if(mystart>(myend_final+mymerge)){  #the new range do not overlap with the current range, "+x" here is a must
      myrow<-matrix(c(mystart_final,myend_final),1,2)  
      mymerge_row<-rbind(mymerge_row,myrow) #store the current range
      mystart_final<-mystart
      myend_final<-myend
      #print("both start and end exceed the range, store the current range and start new")
    }
    if(i==nrow(mytable)){
      myrow<-matrix(c(mystart_final,myend_final),1,2)  
      mymerge_row<-rbind(mymerge_row,myrow) #store the current range
      #print("merging finished")
    }
  }
  
  mymerge_row<-mymerge_row[-1,]
  
  #check the minimum distance between ISs
  mydist<-c()
  for (j in 2:nrow(mymerge_row)){
    myx<-mymerge_row[j,1]-mymerge_row[j-1,2]
    mydist<-c(mydist,myx)
  }
  mymin<-sort(mydist)
  print(mymin[1:10])
  
  mymerge_row<-cbind(as.character(mylist[k]),mymerge_row)
  colnames(mymerge_row)<-c("sseqid","mystart","myend")
  
  mymerge_all<-rbind(mymerge_all,mymerge_row)
}

mymerge_all<-mymerge_all[-1,]

write.table(mymerge_all,file=paste("addIS_merged_size",myISsize,"_freq",myfreq,"_merge",mymerge,".txt",sep=""),quote=F,row.names = F,col.names = T,sep="\t")
  
  
