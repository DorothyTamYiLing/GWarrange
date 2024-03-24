mytable<-read.table("homo_homo_blastout.txt",header=F,sep="\t")
mylist<-unique(c(mytable$V1,mytable$V2))
#both column 1 and column 2 contain all the seqeunce
#length(mylist) #281
#length(unique(mytable$V1)) #281
#length(unique(mytable$V2)) #281

mykeep<-c()
while(length(mylist)>1){
  mycur<-mylist[1]
  mylink<-mytable[which(mytable$V1==mycur),"V2"]
  if(length(mylink)>0){  #if there are some link
    mykeep<-c(mykeep,mycur) #keep the current one
    myrm<-c(mycur,mylink)
    mylist<-mylist[-which(mylist%in%myrm)] #remove the link from the list
    mytable<-mytable[-which(mytable$V1%in%myrm | mytable$V2%in%myrm),] #remove all the rows contain myrm in col1 or col2
  }
  if(length(mylink)==0){ #if no link
    mykeep<-c(mykeep,mycur) #keep the current one
    mylist<-mylist[-which(mylist%in%mycur)]
    mytable<-mytable[-which(mytable$V1%in%mycur | mytable$V2%in%mycur),] #remove 
  }
  if(length(mylist)==1){ #when only one object left
    mykeep<-c(mykeep,mycur) #keep current object
  }
  
}
mykeep<-unique(mykeep)
write.table(mykeep,file="homo_deduplist.txt",sep="\t",row.names=F,col.names = F,quote=F)

