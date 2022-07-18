library("optparse")

option_list = list(
  make_option("--k.len", type="character", default=NULL, 
              help="kmer to plot", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); 


#headers from the blast output
#query, subject, identity, alig_len, mismatches, gap, qstart, qend, sStart, sEnd, evalue, bitscore

###blast output quality check
#all kmer that should show at least one blast hit

#(R)
#load in the blast output file 
mytable<-read.table("myout.txt", header=F)
colnames(mytable)<-c("query","subject","identity","alig_len","mismatches","gap","qstart","qend","sStart","sEnd","evalue","bitscore")

#load in the flank start and end coordinates of the sig kmers
myflk_coor<-read.delim("flank_coor.txt",header=F,sep="_")
colnames(myflk_coor)<-c("kmer","leftflankend","rightflankstart","kmer_len")

#load in phenotype file
myphenofile<-read.table("phenotypes.tsv",header=F)

mykmer<-as.character(unique(mytable$query))  #get the list of kmers with blast output
mygen<-as.character(unique(myphenofile$V1))  #get the list of the genomes from the pheno file

#set the output for the rows of kmers with deletions
del_k<-c()

#set the output for the rows with multiple blast hits in flank 
multi_hit_k<-c()

#set the output for the rows with incomplete flank alignment
alignlen_issue_k<-c()

#set the output for the rows with SNPs and gaps, based on the blast output "mismatches" and "gap" columns 
SNPgap_k<-c()

#create matrix "myprocess" for storing kmers blast hits that are present in all genomes with both flanks; the flanks are also fully aligned with no SNPs nor gaps, and the flanks show unique blast hit in each genome
myprocess<-matrix(0,0,ncol(mytable))
colnames(myprocess)<-c("query","subject","identity","alig_len","mismatches","gap","qstart","qend","sStart","sEnd","evalue","bitscore")

#looping through kmers and genomes
for (i in 1:length(mykmer)){

count=0
myk_row<-mytable[which(mytable$query==mykmer[i]),] #rows referring to the kmer
myflk_coor_k<-unlist(c(1,myflk_coor[which(as.character(myflk_coor$kmer)==as.character(mykmer[i])),2:4]))   #extract the flank start and end coordinate of the kmer from "myflk_coor" file

for (j in 1:length(mygen)){

mysub<-mytable[which(mytable$query==as.character(mykmer[i]) & mytable$subject==as.character(mygen[j])),]
myk_blastcoor<-c(mysub$qstart,mysub$qend)  #extract the blast hit coordinates 
mysnpgap<-c(mysub$mismatches,mysub$gap)

if (nrow(mysub)<=1 | length(unique(myk_blastcoor))<4){  #output the kmers with deleted flank
del_k<-c(del_k,mykmer[i])
count=count+1
}
if (nrow(mysub)>2 | (nrow(mysub)==2 & all(mysub[1,c("qstart","qend")]==mysub[2,c("qstart","qend")]))){   #output the kmers with multi-hit flank
multi_hit_k<-c(multi_hit_k,mykmer[i])
count=count+1
}
if (length(unique(myk_blastcoor))==4 & (any(is.element(myk_blastcoor,myflk_coor_k)!=T | any(is.element(myflk_coor_k,myk_blastcoor)!=T)))){  #output the kmers with incomplete flank blast match
 alignlen_issue_k<-c(alignlen_issue_k,mykmer[i])
 count=count+1
} 
if (any(mysnpgap!=0)) {  #output the kmers with SNP or gaps in flanks (with potential for user to pick the number of SNP and gaps allowed in a flank blast match)
 SNPgap_k<-c(SNPgap_k,mykmer[i])
 count=count+1
} 
 
 } #closing for loop for genomes
 
 if (count==0){
 myprocess<-rbind(myprocess,myk_row)
 }
 
 } #closing for loop for kmers

#the myprocess table should refere to kmers that are present in all genomes with both flanks; the flanks are also fully aligned with no SNPs nor gaps, and the flanks show unique blast hit in each genome 
write.table(myprocess,file="rows_for_process.txt",quote=F,row.names = F,col.names = T,sep="\t")


#output the rows of kmers with quality issues

if (length(unique(del_k))>0){
my_del_k<-mytable[which(mytable$query%in%unique(del_k)),]
write.table(my_del_k,file="kmer_with_deletion.txt",quote=F,row.names = F,col.names = T,sep="\t")
}

if (length(unique(multi_hit_k))>0){
my_multi_hit_k<-mytable[which(mytable$query%in%unique(multi_hit_k)),]
write.table(my_multi_hit_k,file="kmer_with_multi_hits.txt",quote=F,row.names = F,col.names = T,sep="\t")
}

if (length(unique(alignlen_issue_k))>0){
my_alignlen_issue_k<-mytable[which(mytable$query%in%unique(alignlen_issue_k)),]
write.table(my_alignlen_issue_k,file="kmer_with_alignlen_issue.txt",,quote=F,row.names = F,col.names = T,sep="\t")
}

if (length(unique(SNPgap_k))>0){
my_SNPgap_k<-mytable[which(mytable$query%in%unique(SNPgap_k)),]
write.table(my_SNPgap_k,file="kmer_with_SNPgap.txt",quote=F,row.names = F,col.names = T,sep="\t")
}  


###R function determining the flank behaviour 

define_flk_behave<-function(StartL,EndL,StartR,EndR){
	StartL<-as.numeric(as.character(StartL))
	EndL<-as.numeric(as.character(EndL))
	StartR<-as.numeric(as.character(StartR))
	EndR<-as.numeric(as.character(EndR))
	
	#test intact kmer, forward k
     if ((StartL < EndL) & (EndL < StartR) & (StartR < EndR) & ((StartR-EndL) < 2500)){
     mybehave<-"intact_k"
     flk_dist<-(StartR-EndL)
     } else if 
     
     #test intact kmer, reverse k
     ((EndR < StartR) & (StartR < EndL) & (EndL < StartL) & ((EndL-StartR) < 2500)){
     mybehave<-"intact_k"
     flk_dist<-(EndL-StartR)
     } else if
     
     #test flank sequence move away from each other, forward kmer
     ((StartL < EndL) & (EndL < StartR) & (StartR < EndR) & ((StartR-EndL) > 2500)){
     mybehave<-"mv_away"
     flk_dist<-(StartR-EndL)
     } else if
     
     #test flank sequence move away from each other, reverse kmer
     ((StartL > EndL) & (EndL > StartR) & (StartR > EndR) & ((EndL-StartR) > 2500)){
     mybehave<-"mv_away"
     flk_dist<-(EndL-StartR)
     } else if
     
     #test Left flank and right flank swap position, forward kmer
     ((StartR < EndR) & (EndR < StartL) & (StartL < EndL)){
     mybehave<-"swp_flk"
     flk_dist<-(StartL-EndR)
     } else if
     
     #test Left flank and right flank swap position, reverse kmer
     ((EndL < StartL) & (StartL < EndR) & (EndR < StartR)){
     mybehave<-"swp_flk"
     flk_dist<-(EndR-StartL)
     } else if
     
     #test for the presence of inversion
     ((StartL < EndL) & (EndL < EndR) & (EndR < StartR)){
     mybehave<-"mv&flp"
     flk_dist<-(EndR-EndL)
     } else if

     ((EndL < StartL) & (StartL < StartR) & (StartR < EndR)){
     mybehave<-"mv&flp"
     flk_dist<-(StartR-StartL)
     } else if
     
     ((StartR < EndR) & (EndR < EndL) & (EndL < StartL)){
     mybehave<-"mv&flp"
     flk_dist<-(EndL-EndR)
     } else if
     
     ((EndR < StartR) & (StartR < StartL) & (StartL < EndL)){
     mybehave<-"mv&flp"
     flk_dist<-(StartL-StartR)
     } else {

     #output undefined behaviour of the flanking sequences
     mybehave<-"undefined_behave"
     flk_dist<-"NA"
     }
     
     return(c(mybehave,flk_dist))     
}


###determining StartL, EndL, StartR, EndR and flank behaviour for each kmer and genome combination

#need to define k_len as argument
#k_len=200

#creating a new table storing the StartL, EndL, StartR, EndR for each kmer and genome combination blast result
mystartendLR<-matrix(0,0,8)
colnames(mystartendLR)<-c("kmer","genome","StartL","EndL","StartR","EndR","mybehave","flk_dist")

#get the list of kmers with both flanks present
mykmer<-unique(myprocess$query) 

#get the list of genomes in myprocess table
mygen<-unique(myprocess$subject) 

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
myrightflank<-mysub[which(mysub$qstart==opt$k.len | mysub$qend==opt$k.len),]   
if(myrightflank$qstart==opt$k.len){ 
StartR<-myrightflank$sEnd
EndR<-myrightflank$sStart
} 
if(myrightflank$qend==opt$k.len){
StartR<-myrightflank$sStart
EndR<-myrightflank$sEnd
}

#call the function for define flank behaviour, output are mybehave and flk_dist
mybehave=flk_dist=NA #clear the variable
myreturn<-define_flk_behave(StartL,EndL,StartR,EndR) 
mybehave<-myreturn[1]
flk_dist<-myreturn[2]

myoutrow<-as.matrix(c(as.character(mykmer[i]),as.character(mygen[j]),StartL,EndL,StartR,EndR,mybehave,flk_dist),1,8)
myoutrow<-t(myoutrow)
mystartendLR<-rbind(mystartendLR,myoutrow)

} #closing loop through each genome

} #closing loop through each kmer

mystartendLR<-as.data.frame(mystartendLR)


#put all the rows of kmers with at least one undefined behaviour into a table for output
myundefine_out<-matrix(0,0,8)
colnames(myundefine_out)<-colnames(mystartendLR)

myk_undefine<-mystartendLR[which(mystartendLR$mybehave=="undefined_behave"),"kmer"]
myundefine_out<-mystartendLR[which(mystartendLR$kmer%in%myk_undefine),]

write.table(myundefine_out,file="myundefine_k.txt",quote=F,row.names = F,col.names = T,sep="\t")

#select the rows with kmers with no undefined behaviour for merging with phenotype
'%!in%' <- function(x,y)!('%in%'(x,y)) #creating the function
mystartendLR_out<-mystartendLR[which(mystartendLR$kmer%!in%myk_undefine),]

#merge in the phenotype information
myflk_behave_pheno<-merge(myphenofile,mystartendLR_out,by.x="V1",by.y="genome")
myflk_behave_pheno<-myflk_behave_pheno[order(myflk_behave_pheno$kmer, decreasing=T),]
colnames(myflk_behave_pheno)[1]<-"genome"
colnames(myflk_behave_pheno)[2]<-"case_control"

write.table(myflk_behave_pheno,file="myflk_behave_pheno.txt",quote=F,row.names = F,col.names = T,sep="\t")


###making flank behaviour summary table for each kmer across all genomes

myflk_behave_pheno$StartL<-as.numeric(as.character(myflk_behave_pheno$StartL))
myflk_behave_pheno$EndL<-as.numeric(as.character(myflk_behave_pheno$EndL))
myflk_behave_pheno$StartR<-as.numeric(as.character(myflk_behave_pheno$StartR))
myflk_behave_pheno$EndR<-as.numeric(as.character(myflk_behave_pheno$EndR))
myflk_behave_pheno$flk_dist<-as.numeric(as.character(myflk_behave_pheno$flk_dist))


#make the final output
myall_out<-matrix(0,1,13)
colnames(myall_out)<-c("kmer","event_sum","flk_behaviour","case_assos","case_assos_prop","ctrl_assos","ctrl_assos_prop","case_assos_gp_Lflk_sumstat","case_assos_gp_Rflk_sumstat","ctrl_assos_gp_Lflk_sumstat","ctrl_assos_gp_Rflk_sumstat","case_assos_gp_flkdis_sumstat","ctrl_assos_gp_flkdis_sumstat")

#extract the unique kmer
myk4plot<-unique(myflk_behave_pheno$kmer)

for (j in 1:length(myk4plot)){ #open bracket for looping through each kmer
  #for (j in 1:10){ #open bracket for looping through each kmer
  
  mykmer<-myk4plot[j] 
  
  #select the rows referring to the kmer
  mytable<-myflk_behave_pheno[which(myflk_behave_pheno$kmer==mykmer),]
  
  #making the output matrix
  myout<-matrix(0,1,13)
  colnames(myout)<-c("kmer","event_sum","flk_behaviour","case_assos","case_assos_prop","ctrl_assos","ctrl_assos_prop","case_assos_gp_Lflk_sumstat","case_assos_gp_Rflk_sumstat","ctrl_assos_gp_Lflk_sumstat","ctrl_assos_gp_Rflk_sumstat","case_assos_gp_flkdis_sumstat","ctrl_assos_gp_flkdis_sumstat")
  myout<-as.data.frame(myout)
  myout$kmer<-mykmer
  
  #get the total number of cases and controls
  ctrl_count<-length(which(mytable$case_control=="0"))
  case_count<-length(which(mytable$case_control=="1"))
  
  #count the proportion of cases and controls in mybehave column
  mytable$mybehave<-as.character(mytable$mybehave) 
  mycat<-unique(mytable$mybehave)
  
  myout$event_sum<-paste(mycat,collapse=":") #fill in the table
  
  mysum_str<-""  #pasting different behaviours into one string
  
  for (i in 1:length(mycat)){ #looping through each behaviour
    mygp<-as.character(mycat[i])  #extract the behave group name
    
    #count number and proportion of gp in ctrl genomes
    mygp_ctrl<-length(which(mytable$mybehave==mygp & mytable$case_control=="0"))
    mygp_ctrl_prop<-round(mygp_ctrl/ctrl_count,2)
    myctrl_str<-paste(mygp_ctrl,"/",ctrl_count,"(",mygp_ctrl_prop,")",sep="")
    
    #count number and proportion of gp in case genomes
    mygp_case<-length(which(mytable$mybehave==mygp & mytable$case_control=="1"))
    mygp_case_prop<-round(mygp_case/case_count,2)
    mycase_str<-paste(mygp_case,"/",case_count,"(",mygp_case_prop,")",sep="")
    
    mysum<-paste(mygp,mycase_str,myctrl_str,sep=":") #make summary string for each behaviour
    
    mysum_str<-paste(mysum_str,mysum,sep=" ") #pasting different behaviours into one string
    
    #define the case and control associated flank behaviour 
    
    #if this behaviour is associated with ctrl 
    if(mygp_ctrl_prop>0.6 & mygp_case_prop<0.4){   
      myctrl_assos_gp<-mygp  #define this variable for making coordinate summary statistics
      myout$ctrl_assos<-mygp
      myout$ctrl_assos_prop<-mygp_ctrl_prop
      
    #get the coordinate summary statistics (based on StartL and StartR) of the myctrl_assos_gp in ctrl genomes, no information on the flank direction
	
	myStartLstat<-paste(round(summary(mytable[which(mytable$mybehave==myctrl_assos_gp & mytable$case_control=="0"),"StartL"]),0),collapse=" ")
	myStartLSD<-round(sd(mytable[which(mytable$mybehave==myctrl_assos_gp & mytable$case_control=="0"),"StartL"]),0)
	myctrl_Lflk_sumstat_str<-paste(myStartLstat,myStartLSD,sep=" ")

	myStartRstat<-paste(round(summary(mytable[which(mytable$mybehave==myctrl_assos_gp & mytable$case_control=="0"),"StartR"]),0),collapse=" ")
	myStartRSD<-round(sd(mytable[which(mytable$mybehave==myctrl_assos_gp & mytable$case_control=="0"),"StartR"]),0)
	myctrl_Rflk_sumstat_str<-paste(myStartRstat,myStartRSD,sep=" ")

	myctrl_flkdiststat_str<-paste(round(summary(mytable[which(mytable$mybehave==myctrl_assos_gp & mytable$case_control=="0"),"flk_dist"]),0),collapse=" ")

	myout$ctrl_assos_gp_Lflk_sumstat<-myctrl_Lflk_sumstat_str  #fill in the table
	myout$ctrl_assos_gp_Rflk_sumstat<-myctrl_Rflk_sumstat_str  #fill in the table
	myout$ctrl_assos_gp_flkdis_sumstat<-myctrl_flkdiststat_str  #fill in the table

    } 
    
    #if this behaviour is associated with case 
    
    if(mygp_case_prop>0.6 & mygp_ctrl_prop<0.4){ 
      mycase_assos_gp<-mygp    #define this variable for making coordinate summary statistics
      myout$case_assos<-mygp
      myout$case_assos_prop<-mygp_case_prop
      
      #clear the variables
	myStartLstat=myStartLSD=myStartRstat=myStartRSD=myflkdiststat=NA

	#get the coordinate summary statistics (based on StartL and StartR) of the mycase_assos_gp in case genomes, no information on the flank direction

	myStartLstat<-paste(round(summary(mytable[which(mytable$mybehave==mycase_assos_gp & mytable$case_control=="1"),"StartL"]),0),collapse=" ")
	myStartLSD<-round(sd(mytable[which(mytable$mybehave==mycase_assos_gp & mytable$case_control=="1"),"StartL"]),0)
	mycase_Lflk_sumstat_str<-paste(myStartLstat,myStartLSD,sep=" ")

	myStartRstat<-paste(round(summary(mytable[which(mytable$mybehave==mycase_assos_gp & mytable$case_control=="1"),"StartR"]),0),collapse=" ")
	myStartRSD<-round(sd(mytable[which(mytable$mybehave==mycase_assos_gp & mytable$case_control=="1"),"StartR"]),0)
	mycase_Rflk_sumstat_str<-paste(myStartRstat,myStartRSD,sep=" ")

	mycase_flkdiststat_str<-paste(round(summary(mytable[which(mytable$mybehave==mycase_assos_gp & mytable$case_control=="1"),"flk_dist"]),0),collapse=" ")

	myout$case_assos_gp_Lflk_sumstat<-mycase_Lflk_sumstat_str  #fill in the table
	myout$case_assos_gp_Rflk_sumstat<-mycase_Rflk_sumstat_str  #fill in the table
	myout$case_assos_gp_flkdis_sumstat<-mycase_flkdiststat_str  #fill in the table
  
    } 
    
    } #close bracket for looping through each behaviour
    
    myout$flk_behaviour<-mysum_str   #fill in the table with the behaviour summary
    
  myall_out<-rbind(myall_out,myout)
  
} #close bracket for looping through each kmer

write.table(myall_out,file="myall_out.txt",quote=F,row.names = F,col.names = T,sep="\t")










