library("optparse")

option_list = list(
  make_option("--dist", type="character", default=1000,
              help="adjacent homo seq with >=xxxbp distance is defined as separate homo cluster", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


mydist<-as.numeric(as.character(opt$dist))

#mydist<-1000

mytable<-read.table("homo_blastdedup_out.txt",header=F)

colnames(mytable)<-c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
mytable$start1<-NA
mytable$end1<-NA

myfwdhit<-which(mytable$sstart<mytable$send)
mytable$start1[myfwdhit]<-mytable$sstart[myfwdhit]
mytable$end1[myfwdhit]<-mytable$send[myfwdhit]

myrevhit<-which(mytable$sstart>mytable$send)
mytable$start1[myrevhit]<-mytable$send[myrevhit]
mytable$end1[myrevhit]<-mytable$sstart[myrevhit]

#need to order the table by start1
mytable<-mytable[order(mytable$start1,decreasing=F),]
distance<-c()
for (i in 1:nrow(mytable)-1){
  genedist<-mytable$start1[i+1]-mytable$end1[i]
  distance<-c(distance,genedist)
}

#set big value in last row 
distance<-c(distance,"1000000")

mytable<-cbind(mytable,distance)

mytable$distance<-as.numeric(as.character(mytable$distance))
mytable$start1<-as.numeric(as.character(mytable$start1))
mytable$end1<-as.numeric(as.character(mytable$end1))

#detect homo cluster size
mytable$cluster_id<-NA
mytable$cluster_size<-NA
cluster_id<-1
count<-0
for (i in 1:nrow(mytable)){
  if(count==0 & mytable$distance[i]<mydist){ #start of a homo
    myhomostart<-i
    count=count+1
  }
  if(count!=0 & mytable$distance[i]>=mydist){ #end of homo
    myhomoend<-i
    if((myhomoend-myhomostart)>1){ #if there are >=2 homo in the cluster
    mysize<-mytable$end1[myhomoend]-mytable$start1[myhomostart]
    mytable$cluster_id[myhomostart:myhomoend]<-cluster_id
    mytable$cluster_size[myhomostart:myhomoend]<-mysize
    cluster_id<-cluster_id+1
    }
    #reset everything
    count<-0
    myhomostart<-0
    myhomoend<-0
  }
  }

write.table(mytable,file="homo_cluster.txt",quote=F,col.names=T,row.names=F,sep="\t")

mylargest<-max(mytable$cluster_size[!is.na(mytable$cluster_size)])
print(paste("size of largest repeat loci cluster : ",mylargest),sep="")

mycount<-as.data.frame(table(mytable$qseqid))
write.table(mycount,file="homo_occurence.txt",quote=F,col.names=F,row.names=F,sep="\t")
