#Usage: Rscript merge_IS.R --input C505_IS481_IS1663_IS110.txt --extend 100 --merge 200 


library("optparse")

option_list = list(
  make_option("--input", type="character", default=NULL, 
              help="input_file", metavar="character"),
  make_option("--extend", type="character", default=NULL, 
              help="number of bp to extend from the start and end of each IS blast hit", metavar="character"),
  make_option("--merge", type="character", default=NULL, 
              help="maximum number of bp between two adjacent IS for them to be merged into one, this number can be the length of kmer used in GWAS", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); 

myextend<-as.numeric(opt$extend)
mymerge<-as.numeric(opt$merge)

myinput<-read.delim(opt$input,header=F,sep="\t")

myncol<-ncol(myinput)

myinput[ , myncol+1] <- NA  #add the new start column
myinput[ , myncol+2] <- NA  #add the new end column

colnames(myinput)<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","mystart","myend")  

#make qstart to be the smaller value and qend to be the larger
myinput[,"mystart"]<-myinput[,"sstart"]  #13 col always take the smaller value
myinput[,"myend"]<-myinput[,"send"] #14 col always take the bigger value

myreverse<-which(myinput[,"sstart"]>myinput[,"send"])
myinput[myreverse,"mystart"]<-myinput[myreverse,"send"]  #13 col always take the smaller value
myinput[myreverse,"myend"]<-myinput[myreverse,"sstart"] #14 col always take the bigger value

#then add 00bp to each side of the IS coordinate to make sure that the whole IS is covered
myinput[,"mystart"]<-myinput[,"mystart"]-myextend
myinput[,"myend"]<-myinput[,"myend"]+myextend

mymerge_all<-matrix(0,1,3)
colnames(mymerge_all)<-c("sseqid","mystart","myend")

#loop by samples
mylist<-unique(myinput$sseqid)

for(k in 1:length(mylist)){
  print(mylist[k])
  mytable<-myinput[which(myinput$sseqid==mylist[k]),]
  
  #mytable<-myinput[which(myinput$sseqid=="F578"),]
  
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

write.table(mymerge_all,file="all_mergedIS.txt",quote=F,row.names = F,col.names = T,sep="\t")

