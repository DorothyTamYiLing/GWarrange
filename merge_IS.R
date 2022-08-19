#usage:Rscript merge_IS.R -i myblastout.txt

#IS replacement
library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="file containing the list of samples to make summary of their local ancestry fragments blast output, in a .txt file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); 

myinput<-read.table(opt$input,header=F)
#myinput<-read.table("/Users/dorothytam/Desktop/myblastout.txt",header=F)

myinput[ , 13] <- NA  #add the new start column
myinput[ , 14] <- NA  #add the new end column
myinput[ , 15] <- NA  #add the merge info column

colnames(myinput)<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","mystart","myend","merge")

myallout<-matrix(0,1,15)

colnames(myallout)<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","mystart","myend","merge")

mylist<-unique(myinput$sseqid)

for(i in 1:length(mylist)){
  print(mylist[i])
  mytable<-myinput[which(myinput$sseqid==mylist[i]),]
  
  #make qstart to be the smaller value and qend ro be the larger
  for (m in 1:nrow(mytable)){
    if(mytable[m,9]<mytable[m,10]){
      mytable[m,13]<-mytable[m,9]  #13 col always take the smaller value
      mytable[m,14]<-mytable[m,10] #14 col always take the bigger value
    }else{
      mytable[m,13]<-mytable[m,10]
      mytable[m,14]<-mytable[m,9]
    }
  }
  
  #find all the rows that are involved in overlapping IS elements
  all_dist<-c()
  mytable<-mytable[order(mytable$mystart),]
  for(j in 2:nrow(mytable)){
    mystart<-mytable[j,13]
    myend<-mytable[j-1,14]
    mydist<-mystart-myend
    all_dist<-c(all_dist,mydist)  #save all the distance values that are negative
  }
  
  myneg_row<-which(all_dist<0 | all_dist==0)
  
  
  #mytable_ori<- mytable # for checking

  #now process those row pairs that involved in each negative distance 
  #(i.e.index of the negative distance value & that +1)
  #In order to take into account of the possibility of 3 consecutive IS, 
  #I transfer the end coordinate from the second row in the pair onto the first row, update the "hit number" column,
  #then delete the second row
  
  for (i in 1:length(myneg_row)){
    myrow<-myneg_row[i]
    mytable[myrow,14]<-mytable[myrow+1,14]  #transfer the end coordinate from second row to first row
    #mytable[myrow,4]<-paste(mytable[myrow,4],"and",mytable[myrow+1,4], sep="") #update the "hit number" col in the fist row
    #mytable[myrow,13]<-paste(mytable[myrow,13],mytable[myrow+1,13], sep="&") #update the "orientation" col in the fist row
    mytable[myrow,15]<-"merged_IS" #remove unnecessary info
    mytable<-mytable[-(myrow+1),] #remove the second row in the pair
    myneg_row<-myneg_row-1 #update the row index since one row is deleted
  }
  
myallout<-rbind(myallout,mytable)
}

myallout<-myallout[-1,]
myallout<-myallout[,c(2,13,14,15)]

write.table(myallout,file="myblastout_mergedIS.txt",quote=F,col.names = T,row.names = F,sep="\t")

