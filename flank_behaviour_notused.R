#determinin StartL, EndL, StartR, EndR for each kmer and genome combination
mykmer<-unique(myprocess$query) #get the list of kmers with both flanks present

k_len=200

#creating a new table storing the StartL, EndL, StartR, EndR for each kmer and genome combination blast result
mystartendLR<-matrix(0,0,6)
colnames(mystartendLR)<-c("kmer","genome","StartL","EndL","StartR","EndR")

#loop through each kmer
for (i in 1:length(mykmer)){
#loop through each genome, and 
for (j in 1:length(mygen)){
#extract the rows referring to the kmer and genome combination
mysub<-myprocess[which(myprocess$query==as.character(mykmer[i]) & myprocess$subject==as.character(mygen[j])),]
#determining left flank, hence StartL and EndL
myleftflank<-mysub[which(mysub$qstart==1 | mysub$qend==1),]  
if(myleftflank$qstart==1){ 
StartL<-myleftflank$sStart
EndL<-myleftflank$sEnd
} 
if(myleftflank$qend==1){
StartL<-myleftflank$sEnd
EndL<-myleftflank$sStart
}
#determining right flank, hence StartR and EndR
myrightflank<-mysub[which(mysub$qstart==k_len | mysub$qend==k_len),]   
if(myrightflank$qstart==k_len){ 
StartR<-myrightflank$sEnd
EndR<-myrightflank$sStart
} 
if(myrightflank$qend==k_len){
StartR<-myrightflank$sStart
EndR<-myrightflank$sEnd
}
myoutrow<-as.matrix(c(as.character(mykmer[i]),mygen[j],StartL,EndL,StartR,EndR),1,6)
myoutrow<-t(myoutrow)
mystartendLR<-rbind(mystartendLR,myoutrow)
} #closing loop through each genome
} #closing loop through each kmer

mystartendLR<-as.data.frame(mystartendLR)

####################################################

#determining the flank behaviour and making the all_behaviour_summary.txt file

mystartendLR$StartL<-as.numeric(as.character(mystartendLR$StartL))
mystartendLR$EndL<-as.numeric(as.character(mystartendLR$EndL))
mystartendLR$StartR<-as.numeric(as.character(mystartendLR$StartR))
mystartendLR$EndR<-as.numeric(as.character(mystartendLR$EndR))

mykmer<-unique(mystartendLR$kmer) #get the list of kmers with both flanks present

#set the output for undefined behaviour of the flanking sequences
my_new_case<-matrix(0,0,ncol(mystartendLR))
colnames(my_new_case)<-colnames(mystartendLR)

#set the output matrix for flank_behaviour.txt file
myflk_behave<-matrix(0,0,10)
colnames(myflk_behave)<-c("kmer","genome","flk_behaviour","notes","StartL","EndL","StartR","EndR","flk_dist","myk_orien")

###start of the main loop###

#loop through each kmer
for (i in 1:length(mykmer)){

#get the rows referring to the kmer
myk_row<-mystartendLR[which(mystartendLR$kmer==mykmer[i]),]

#first, determine if the intact kmers are forward or reverse and store in variable "myk_orien", assuming that there will be rows referring to intact kmers

myk_orien=NA  #clear the variable


myintact_fwd<-length(which((myk_row$StartL < myk_row$EndL) & (myk_row$EndL < myk_row$StartR) & (myk_row$StartR < myk_row$EndR) & ((myk_row$StartR-myk_row$EndL) < 2200)))
myintact_rev<-length(which((myk_row$EndR < myk_row$StartR) & (myk_row$StartR < myk_row$EndL) & (myk_row$EndL < myk_row$StartL) & ((myk_row$EndL-myk_row$StartR) < 2200)))
if(myintact_fwd>0 & myintact_rev==0){
myk_orien<-"fwd_k"
} else if
(myintact_fwd==0 & myintact_rev>0){
myk_orien<-"rev_k"
} else {
myk_orien<-"intactk_only"
print("kmer map both forward and reverse in genomes, intact kmers only cases")
#write.table(myk_row,file="both_fwd_rev_intactk.txt",quote=F,row.names = F,col.names = T,sep="\t")
}
 
 #now detecting the flank behaviours, loop through each genome that are in mystartendLR
   
  #looping through each genome
  for (j in 1:length(mygen)){    
     
     mybehave=mynote=StartL=EndL=StartR=EndR=flk_dist=NA #clear the variable
     
     StartL<-myk_row[j,"StartL"]
     EndL<-myk_row[j,"EndL"]
     StartR<-myk_row[j,"StartR"]
     EndR<-myk_row[j,"EndR"]
     
     #test for behaviour 1, intact kmer, forward kmer
     if ((StartL < EndL) & (EndL < StartR) & (StartR < EndR) & ((StartR-EndL) < 2200)){
     mybehave<-"intact_k"
     mynote<-"fwd_k" #store the kmer orientation, the four numbers and the distance between the two flanks
     flk_dist<-(StartR-EndL)
     } else if 
     
     #test for behaviour 1, intact kmer, reverse kmer
     ((EndR < StartR) & (StartR < EndL) & (EndL < StartL) & ((EndL-StartR) < 2200)){
     mybehave<-"intact_k"
     mynote<-"rev_k" #store the kmer orientation, the four numbers and the distance between the two flanks
     flk_dist<-(StartL-EndL)
     } else if
     
     #test for behaviour 2, Flank sequence move away from each other, forward kmer
     ((StartL < EndL) & (EndL < StartR) & (StartR < EndR) & ((StartR-EndL) > 2200)){
     mybehave<-"mv_away"
     mynote<-"fwd_k" #store the kmer orientation, the four numbers and the distance between the two flanks
     flk_dist<-(EndR-StartR)
     } else if
     
     #test for behaviour 2, Flank sequence move away from each other, reverse kmer
     ((StartL > EndL) & (EndL < StartR) & (StartR > EndR) & ((EndL-StartR) > 2200)){
     mybehave<-"mv_away"
     mynote<-"rev_k" #store the kmer orientation, the four numbers and the distance between the two flanks
     flk_dist<-(StartR-EndR)
     } else if
     
     #test for behaviour 3, Left flank and right flank swap position, forward kmer
     ((StartR < EndR) & (EndR < StartL)  & (StartL < EndL)){
     mybehave<-"swp_flk"
     mynote<-"fwd_k"
     flk_dist<-(StartL-EndR)
     } else if
     
     #test for behaviour 3, Left flank and right flank swap position, reverse kmer
     ((EndL < StartL) & (StartL < EndR) & (EndR < StartR)){
     mybehave<-"swp_flk"
     mynote<-"rev_k"
     flk_dist<-(EndR-StartL)
     } else if
     
     #test for behaviour 4 and 5, first look at if the intact k is forward or reverse
     ((StartL < EndL) & (EndL < EndR) & (EndR < StartR) & (myk_orien=="fwd_k")){
     mybehave<-"mv&flp"
     mynote<-"Right_mv&flp"
     flk_dist<-(EndR-EndL)
     } else if

     ((EndL < StartL) & (StartL < StartR) & (StartR < EndR) & (myk_orien=="fwd_k")){
     mybehave<-"mv&flp"
     mynote<-"Left_mv&flp"
     flk_dist<-(StartR-StartL)
     } else if
     
     ((StartR < EndR) & (EndR < EndL) & (EndL < StartL) & (myk_orien=="rev_k")){
     mybehave<-"mv&flp"
     mynote<-"Right_mv&flp"
     flk_dist<-(EndL-EndR)
     } else if
     
     ((EndR < StartR) & (StartR < StartL) & (StartL < EndL) & (myk_orien=="rev_k")){
     mybehave<-"mv&flp"
     mynote<-"Left_mv&flp"
     flk_dist<-(StartL-StartR)
     } else if
     
     ((StartR < EndR) & (EndR < EndL) & (EndL < StartL) & (myk_orien=="fwd_k")){
     mybehave<-"swp&flp"
     mynote<-"Left_swp&flp"
     flk_dist<-(EndL-EndR)
     } else if
     
     ((EndR < StartR) & (StartR < StartL) & (StartL < EndL) & (myk_orien=="fwd_k")){
     mybehave<-"swp&flp"
     mynote<-"Right_swp&flp"
     flk_dist<-(StartL-StartR)
     } else if
     
     ((EndL < StartL) & (StartL < StartR) & (StartR < EndR) & (myk_orien=="rev_k")){
     mybehave<-"swp&flp"
     mynote<-"Right_swp&flp"
     flk_dist<-(StartR-StartL)
     } else if
     
     ((StartL < EndL) & (EndL < EndR) & (EndR < StartR) & (myk_orien=="rev_k")){
     mybehave<-"swp&flp"
     mynote<-"Left_swp&flp"
     flk_dist<-(EndR-EndL)
     } else {

     #output undefined behaviour of the flanking sequences
     my_new_case<-rbind(my_new_case, myk_row[j,])
     }
     
     myoutrow<-c(as.character(mykmer[i]), mygen[j],mybehave,mynote,StartL,EndL,StartR,EndR,flk_dist,myk_orien)
   myflk_behave<-rbind(myflk_behave,myoutrow)  #append the info for output
 
} #closing bracket for the genome for-loop

} #closing brackets for looping through each kmer

###end of the main loop###

#write.table(myflk_behave,file="all_behaviour_summary.txt",quote=F,row.names = F,col.names = T,sep="\t")

#output the kmer rows with undefined flank behaviour
if (nrow(my_new_case)>0){
write.table(my_new_case,file="new_flk_behaviour.txt",quote=F,row.names = F,col.names = T,sep="\t")
}


#load in the phenotype file for merging in the case/control information
myphenofile<-read.delim("case_control.txt", header=F,sep="\t")

#merge in the phenotype information
myflk_behave_pheno<-merge(myphenofile,myflk_behave,by.x="V1",by.y="genome")
myflk_behave_pheno<-myflk_behave_pheno[order(myflk_behave_pheno$kmer, decreasing=T),]

write.table(myflk_behave_pheno,file="myflk_behave_pheno.txt",quote=F,row.names = F,col.names = T,sep="\t")


