library("optparse")
library(plyr)

option_list = list(
  make_option("--k.len", type="character", default=NULL, 
              help="kmer to plot", metavar="character"),
  make_option("--pheno", type="character", default=NULL,
              help="kmer to plot", metavar="character"),
  make_option("--outdir", type="character", default=NULL,
              help="kmer to plot", metavar="character"),
  make_option("--flkdist", type="character", default=NULL,
              help="kmer to plot", metavar="character"),
  make_option("--dedupk", type="character", default=NULL,
              help="kmer to plot", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); 


#headers from the blast output
#query, subject, identity, alig_len, mismatches, gap, qstart, qend, sStart, sEnd, evalue, bitscore

###blast output quality check
#all kmer that should show at least one blast hit

#setwd(opt$outdir)

#(R)
#load in the blast output file 
mytable<-read.table(paste(opt$outdir,"myout.txt",sep="/"), header=F)
colnames(mytable)<-c("query","subject","identity","alig_len","mismatches","gap","qstart","qend","sStart","sEnd","evalue","bitscore")

#load in the flank start and end coordinates of the sig kmers
myflk_coor<-read.delim(paste(opt$outdir,"flank_coor.txt",sep="/"),header=F,sep="_")
colnames(myflk_coor)<-c("kmer","leftflankend","rightflankstart","kmer_len")

#load in phenotype file
myphenofile<-read.table(opt$pheno,header=F)

mykmer<-as.character(unique(mytable$query))  #get the list of kmers with blast output
mygen<-as.character(unique(myphenofile$V1))  #get the list of the genomes from the pheno file

#set the output for the rows of kmers with absence of >5% genomes
abs_gen_k<-c()

#set the output for the rows of kmers with deletions
del_k<-c()

#set the output for the rows with multiple blast hits in flank 
multi_hit_k<-c()

#set the output for the rows with incomplete flank alignment
alignlen_issue_k<-c()

#looping through kmers 
for (i in 1:length(mykmer)){

#print(mykmer[i])
myflk_coor_k<-unlist(c(1,myflk_coor[which(as.character(myflk_coor$kmer)==as.character(mykmer[i])),2:4]))   #extract the flank start and end coordinate of the kmer from "myflk_coor" file
mykrow<-mytable[which(mytable$query==as.character(mykmer[i])),]
#mykrow<-mytable[which(mytable$query=="kmer51"),]


#blast hit must contain >=95% of the genomes used
if(length(unique(mykrow$subject))<(0.95*length(mygen))){
abs_gen_k<-c(abs_gen_k,mykmer[i])
}

#each genome must appear twice
myfreqtable<-data.frame(table(mykrow$subject))
if(any(myfreqtable$Freq<2)){
del_k<-c(del_k,mykmer[i])
}

if(any(myfreqtable$Freq>2)){
multi_hit_k<-c(multi_hit_k,mykmer[i])
}

#kmer blast match coordinates must match , not use this filter anymore, 27/7/2023
#mycoor<-unique(c(mykrow$qstart,mykrow$qend))
#if((any(mycoor%in%myflk_coor_k==F) | any(myflk_coor_k%in%mycoor==F))){
#alignlen_issue_k<-c(alignlen_issue_k,mykmer[i])
#}

#use this filter instead, 27/7/2023
if(any(mykrow$identity<95) | any(mykrow$evalue>0.00005)){
  alignlen_issue_k<-c(alignlen_issue_k,mykmer[i])
}

}#close the for loop

#make summary for kmer quality control
myfilterout<-matrix(0,length(mykmer),4)
rownames(myfilterout)<-mykmer
colnames(myfilterout)<-c("abs_gen_k","del_k","multi_hit_k","alignlen_issue_k")

myfilterout<-as.data.frame(myfilterout)
#myfilterout

#output the rows of kmers with quality issues

if (length(unique(abs_gen_k))>0){
print(paste("abs_gen_k has ",length(unique(abs_gen_k))," kmers",sep=""))
myfilterout[abs_gen_k,"abs_gen_k"]<-"yes"
my_abs_gen_k<-mytable[which(mytable$query%in%unique(abs_gen_k)),]
write.table(my_abs_gen_k,file=paste(opt$outdir,"kmer_with_missinggenomes.txt",sep="/"),quote=F,row.names=F,col.names = T,sep="\t")
}else{
print("no abs_gen_k")
}

if (length(unique(del_k))>0){
print(paste("del_k has ",length(unique(del_k))," kmers",sep=""))
myfilterout[del_k,"del_k"]<-"yes"
my_del_k<-mytable[which(mytable$query%in%unique(del_k)),]
write.table(my_del_k,file=paste(opt$outdir,"kmer_with_missinggenomes.txt",sep="/"),quote=F,row.names=F,col.names = T,sep="\t")
}else{
print("no del_k")
}

if (length(unique(multi_hit_k))>0){
print(paste("multi_hit_k has ",length(unique(multi_hit_k))," kmers",sep=""))
myfilterout[multi_hit_k,"multi_hit_k"]<-"yes"
my_multi_hit_k<-mytable[which(mytable$query%in%unique(multi_hit_k)),]
write.table(my_multi_hit_k,file=paste(opt$outdir,"kmer_with_multi_hits.txt",sep="/"),quote=F,row.names=F,col.names = T,sep="\t")
}else{
print("no multi_hit_k")
}

if (length(unique(alignlen_issue_k))>0){
print(paste("alignlen_issue_k has ",length(unique(alignlen_issue_k))," kmers",sep=""))
myfilterout[alignlen_issue_k,"alignlen_issue_k"]<-"yes"
my_alignlen_issue_k<-mytable[which(mytable$query%in%unique(alignlen_issue_k)),]
write.table(my_alignlen_issue_k,file=paste(opt$outdir,"kmer_with_alignlen_issue.txt",sep="/"),quote=F,row.names=F,col.names = T,sep="\t")
}else{
print("no alignlen_issue_k")
}

write.table(myfilterout,file=paste(opt$outdir,"filterk_out_summary.txt",sep="/"),quote=F,row.names = T,col.names =T,sep="\t")


#output the good kmers for further processing
mybadk<-unique(c(abs_gen_k,del_k,multi_hit_k,alignlen_issue_k))
#mybadk<-unique(c(del_k,multi_hit_k))
mygoodk<-mykmer[which(!is.element(mykmer,mybadk))]

if(length(mygoodk)==0){
print("no good kmers for process")
}else{
print("now process good kmers")

myprocess<-mytable[which(mytable$query%in%mygoodk),]

#the myprocess table should refere to kmers that are present in all genomes with both flanks; the flanks are also fully aligned with no SNPs nor gaps, and the flanks show unique blast hit in each genome 
write.table(myprocess,file=paste(opt$outdir,"rows_for_process.txt",sep="/"),quote=F,row.names = F,col.names = T,sep="\t")

myprocess<-read.table(paste(opt$outdir,"rows_for_process.txt",sep="/"),sep="\t",header=T)

#make a label for kmer:gen combination
myprocess$label<-paste(myprocess$query,myprocess$subject,sep="_")

#extract odd rows
#create a dummy indicator that shows whether a row is even or odd.
row_odd <- seq_len(nrow(myprocess)) %% 2

#then use our dummy to drop all even rows from our data frame
data_row_odd <- myprocess[row_odd == 1, ]   
colnames(data_row_odd)<-c("query_o","subject_o","identity_o","alig_len_o","mismatches_o","gap_o","qstart_o","qend_o","sStart_o","sEnd_o","evalue_o","bitscore_o","label_o")

#create even rows
data_row_even <- myprocess[row_odd == 0, ]   
colnames(data_row_even)<-c("query_e","subject_e","identity_e","alig_len_e","mismatches_e","gap_e","qstart_e","qend_e","sStart_e","sEnd_e","evalue_e","bitscore_e","label_e")

#paste the odd and even rows side by side
mymerge<-merge(data_row_odd,data_row_even,by.x="label_o",by.y="label_e")

#checking rows referring to the same genomes and the same kmer, for myself only
all(mymerge$query_o==mymerge$query_e)
all(mymerge$subject_o==mymerge$subject_e)

#determining StartL, EndL, StartR, EndR for each kmer and genome combination blast result

print("determine StartL, EndL, StartR, EndR")

mymerge$StartL<-0
mymerge$EndL<-0
mymerge$StartR<-0
mymerge$EndR<-0

#################### new ways based on flexible kmers length ###############

for(m in 1:nrow(mymerge)){
myq<-c(mymerge$qstart_o[m],mymerge$qend_o[m],mymerge$qstart_e[m],mymerge$qend_e[m])
  mys<-c(mymerge$sStart_o[m],mymerge$sEnd_o[m],mymerge$sStart_e[m],mymerge$sEnd_e[m])
  #define startL, startR, endL, endR for subject
  mys<-mys[order(myq)] #order subject coor by query coor
  mymerge$StartL[m]<-mys[1]
  mymerge$EndL[m]<-mys[2]
  mymerge$StartR[m]<-mys[3]
  mymerge$EndR[m]<-mys[4]
  #define startL, startR, endL, endR for query
  #myq<-myq[order(myq)]
  #mymerge$startL_q[m]<-myq[1]
  #mymerge$startR_q[m]<-myq[2]
  #mymerge$endL_q[m]<-myq[3]
  #mymerge$endR_q[m]<-myq[4]
}

############################################################################

#################### old ways based on fixed kmers length ##################

#k.len=opt$k.len
#k.len=200

#defining StartL and EndL


#for rows where qstart_o==1
#myqstart_o_1<-which(mymerge$qstart_o==1)
#mymerge[myqstart_o_1,"StartL"]<-mymerge[myqstart_o_1,"sStart_o"] 
#mymerge[myqstart_o_1,"EndL"]<-mymerge[myqstart_o_1,"sEnd_o"] 

#for rows where qend_o==1
#myqend_o_1<-which(mymerge$qend_o==1)
#mymerge[myqend_o_1,"StartL"]<-mymerge[myqend_o_1,"sEnd_o"] 
#mymerge[myqend_o_1,"EndL"]<-mymerge[myqend_o_1,"sStart_o"] 

#for rows where qstart_e==1
#myqstart_e_1<-which(mymerge$qstart_e==1)
#mymerge[myqstart_e_1,"StartL"]<-mymerge[myqstart_e_1,"sStart_e"] 
#mymerge[myqstart_e_1,"EndL"]<-mymerge[myqstart_e_1,"sEnd_e"] 

#for rows where qend_e==1
#myqend_e_1<-which(mymerge$qend_e==1)
#mymerge[myqend_e_1,"StartL"]<-mymerge[myqend_e_1,"sEnd_e"] 
#mymerge[myqend_e_1,"EndL"]<-mymerge[myqend_e_1,"sStart_e"] 


#defining StartR and EndR
#for rows where qstart_o==k.len
#myqstart_o_klen<-which(mymerge$qstart_o==k.len)
#mymerge[myqstart_o_klen,"EndR"]<-mymerge[myqstart_o_klen,"sStart_o"] 
#mymerge[myqstart_o_klen,"StartR"]<-mymerge[myqstart_o_klen,"sEnd_o"] 

#for rows where qend_o==k.len
#myqend_o_klen<-which(mymerge$qend_o==k.len)
#mymerge[myqend_o_klen,"EndR"]<-mymerge[myqend_o_klen,"sEnd_o"] 
#mymerge[myqend_o_klen,"StartR"]<-mymerge[myqend_o_klen,"sStart_o"] 

#for rows where qstart_e==k.len
#myqstart_e_klen<-which(mymerge$qstart_e==k.len)
#mymerge[myqstart_e_klen,"EndR"]<-mymerge[myqstart_e_klen,"sStart_e"] 
#mymerge[myqstart_e_klen,"StartR"]<-mymerge[myqstart_e_klen,"sEnd_e"] 

#for rows where qend_e==k.len
#myqend_e_klen<-which(mymerge$qend_e==k.len)
#mymerge[myqend_e_klen,"StartR"]<-mymerge[myqend_e_klen,"sStart_e"] 
#mymerge[myqend_e_klen,"EndR"]<-mymerge[myqend_e_klen,"sEnd_e"] 

########################################################################

#defining the dist to define split flank, find out the maximum IS replacement block size, then add few thousands bp
#myflkdist=70000 # merge7000_ext500, max IS replacment block size 45118bp
#myflkdist=200000 # merge15000_ext500, max IS replacement block size 159836bp

myflkdist<-as.numeric(as.character(opt$flkdist))


mymerge$mybehave<-0 #creating new columns
mymerge$flk_dist<-0 #creating new columns

mymerge$StartL<-as.numeric(as.character(mymerge$StartL))
mymerge$EndL<-as.numeric(as.character(mymerge$EndL))
mymerge$StartR<-as.numeric(as.character(mymerge$StartR))
mymerge$EndR<-as.numeric(as.character(mymerge$EndR))

print("determine behaviour")

#First, set all rows as "undefined behaviour" and replacing by different condition
mymerge$mybehave<-"undefined_behave"
mymerge$intactk_orien<-"NA"
mymerge$flk_dist<-"NA"

#test intact kmer, forward k
myintactk_fwd_row<-which((mymerge$StartL < mymerge$EndL) & (mymerge$EndL < mymerge$StartR) & (mymerge$StartR < mymerge$EndR) & ((mymerge$StartR-mymerge$EndL) < myflkdist))
mymerge[myintactk_fwd_row,"mybehave"]<-"intact_k"
mymerge[myintactk_fwd_row,"intactk_orien"]<-"intactk_fwd"
mymerge[myintactk_fwd_row,"flk_dist"]<-mymerge[myintactk_fwd_row,"StartR"]-mymerge[myintactk_fwd_row,"EndL"]

#test intact kmer, reverse k
myintactk_rev_row<-which((mymerge$EndR < mymerge$StartR) & (mymerge$StartR < mymerge$EndL) & (mymerge$EndL < mymerge$StartL) & ((mymerge$EndL-mymerge$StartR) < myflkdist))
mymerge[myintactk_rev_row,"mybehave"]<-"intact_k"
mymerge[myintactk_rev_row,"intactk_orien"]<-"intactk_rev"
mymerge[myintactk_rev_row,"flk_dist"]<-mymerge[myintactk_rev_row,"EndL"]-mymerge[myintactk_rev_row,"StartR"]

#test flank sequence move away from each other, forward kmer
mymv_away_fwd_row<-which((mymerge$StartL < mymerge$EndL) & (mymerge$EndL < mymerge$StartR) & (mymerge$StartR < mymerge$EndR) & ((mymerge$StartR-mymerge$EndL) > myflkdist))
mymerge[mymv_away_fwd_row,"mybehave"]<-"mv_away"
mymerge[mymv_away_fwd_row,"flk_dist"]<-mymerge[mymv_away_fwd_row,"StartR"]-mymerge[mymv_away_fwd_row,"EndL"]

#test flank sequence move away from each other, reverse kmer
mymv_away_rev_row<-which((mymerge$StartL > mymerge$EndL) & (mymerge$EndL > mymerge$StartR) & (mymerge$StartR > mymerge$EndR) & ((mymerge$EndL-mymerge$StartR) > myflkdist))
mymerge[mymv_away_rev_row,"mybehave"]<-"mv_away"
mymerge[mymv_away_rev_row,"flk_dist"]<-mymerge[mymv_away_rev_row,"EndL"]-mymerge[mymv_away_rev_row,"StartR"]

#test Left flank and right flank swap position, forward kmer
myswp_flk_fwd_row<-which((mymerge$StartR < mymerge$EndR) & (mymerge$EndR < mymerge$StartL) & (mymerge$StartL < mymerge$EndL))
mymerge[myswp_flk_fwd_row,"mybehave"]<-"swp_flk"
mymerge[myswp_flk_fwd_row,"flk_dist"]<-mymerge[myswp_flk_fwd_row,"StartL"]-mymerge[myswp_flk_fwd_row,"EndR"]

#test Left flank and right flank swap position, reverse kmer
myswp_flk_rev_row<-which((mymerge$EndL < mymerge$StartL) & (mymerge$StartL < mymerge$EndR) & (mymerge$EndR < mymerge$StartR))
mymerge[myswp_flk_rev_row,"mybehave"]<-"swp_flk"
mymerge[myswp_flk_rev_row,"flk_dist"]<-mymerge[myswp_flk_rev_row,"EndR"]-mymerge[myswp_flk_rev_row,"StartL"]

#test for the presence of inversion
my_mvandflp_1_row<-which((mymerge$StartL < mymerge$EndL) & (mymerge$EndL < mymerge$EndR) & (mymerge$EndR < mymerge$StartR))
mymerge[my_mvandflp_1_row,"mybehave"]<-"mv&flp"
mymerge[my_mvandflp_1_row,"flk_dist"]<-mymerge[my_mvandflp_1_row,"EndR"]-mymerge[my_mvandflp_1_row,"EndL"]
    
my_mvandflp_2_row<-which((mymerge$EndL < mymerge$StartL) & (mymerge$StartL < mymerge$StartR) & (mymerge$StartR < mymerge$EndR))
mymerge[my_mvandflp_2_row,"mybehave"]<-"mv&flp"
mymerge[my_mvandflp_2_row,"flk_dist"]<-mymerge[my_mvandflp_2_row,"StartR"]-mymerge[my_mvandflp_2_row,"StartL"]
 
my_mvandflp_3_row<-which((mymerge$StartR < mymerge$EndR) & (mymerge$EndR < mymerge$EndL) & (mymerge$EndL < mymerge$StartL))
mymerge[my_mvandflp_3_row,"mybehave"]<-"mv&flp"    
mymerge[my_mvandflp_3_row,"flk_dist"]<-mymerge[my_mvandflp_3_row,"EndL"]-mymerge[my_mvandflp_3_row,"EndR"]

my_mvandflp_4_row<-which((mymerge$EndR < mymerge$StartR) & (mymerge$StartR < mymerge$StartL) & (mymerge$StartL < mymerge$EndL))
mymerge[my_mvandflp_4_row,"mybehave"]<-"mv&flp"    
mymerge[my_mvandflp_4_row,"flk_dist"]<-mymerge[my_mvandflp_4_row,"StartL"]-mymerge[my_mvandflp_4_row,"StartR"]

#selecting columns StartL, EndL, StartR, EndR for each kmer and genome combination blast result
mymerge<-mymerge[,c("query_o","subject_o","StartL","EndL","StartR","EndR","mybehave","intactk_orien","flk_dist")]
colnames(mymerge)<-c("kmer","genome","StartL","EndL","StartR","EndR","mybehave","intactk_orien","flk_dist")

mystartendLR<-as.data.frame((mymerge))

#myphenofile<-read.table("../prn_status_pheno.txt",header=F)

#put all the rows of kmers with at least one undefined behaviour into a table for output
myk_undefine<-mystartendLR[which(mystartendLR$mybehave=="undefined_behave"),"kmer"]

if (length(unique(myk_undefine))>0){
myundefine_out<-matrix(0,0,8)
colnames(myundefine_out)<-colnames(mystartendLR)
myundefine_out<-mystartendLR[which(mystartendLR$kmer%in%myk_undefine),]
write.table(myundefine_out,file="myundefine_k.txt",quote=F,row.names = F,col.names = T,sep="\t")
}

#select the rows with kmers with no undefined behaviour for merging with phenotype
'%!in%' <- function(x,y)!('%in%'(x,y)) #creating the function
mystartendLR_out<-mystartendLR[which(mystartendLR$kmer%!in%myk_undefine),]

#merge in the phenotype information
myflk_behave_pheno<-merge(myphenofile,mystartendLR_out,by.x="V1",by.y="genome")
myflk_behave_pheno<-myflk_behave_pheno[order(myflk_behave_pheno$kmer, decreasing=T),]
colnames(myflk_behave_pheno)[1]<-"genome"
colnames(myflk_behave_pheno)[2]<-"case_control"

write.table(myflk_behave_pheno,file=paste(opt$outdir,"myflk_behave_pheno.txt",sep="/"),quote=F,row.names = F,col.names = T,sep="\t")


###making flank behaviour summary table for each kmer across all genomes
#First, make the summary only for kmers with flank behaviour other than "intact_k"

print("make summary for split k")

#read in myflk_behave_pheno.txt
#myall_behave<-read.table("myflk_behave_pheno.txt",sep="\t",header=T)
myall_behave<-myflk_behave_pheno

'%!in%' <- function(x,y)!('%in%'(x,y)) #creating the function

#list of kmers with flank behaviour other than "intact_k"
myk_otherbebave<-unique(myall_behave[which(myall_behave$mybehave!="intact_k"),"kmer"])

if(length(myk_otherbebave)>0){

#rows of kmers with flank behaviour other than "intact_k"
myflk_behave_pheno<-myall_behave[which(myall_behave$kmer%in%myk_otherbebave),]

#then run the lines for making summary for kmers with flank behaviour other than "intact_k"
myflk_behave_pheno$StartL<-as.numeric(as.character(myflk_behave_pheno$StartL))
myflk_behave_pheno$EndL<-as.numeric(as.character(myflk_behave_pheno$EndL))
myflk_behave_pheno$StartR<-as.numeric(as.character(myflk_behave_pheno$StartR))
myflk_behave_pheno$EndR<-as.numeric(as.character(myflk_behave_pheno$EndR))
myflk_behave_pheno$flk_dist<-as.numeric(as.character(myflk_behave_pheno$flk_dist))


#make the final output
myall_out<-matrix(0,1,9)
colnames(myall_out)<-c("kmer","event_sum","flk_behaviour","my0_intactk_sum","my1_intactk_sum","otherk","my0_otherk_sum","my1_otherk_sum","event")

#make the final brief output for signifiant behaviours only
myshort_allout<-matrix(0,1,8)
colnames(myshort_allout)<-c("kmer","intactk_mygp_ctrl_prop","intactk_mygp_case_prop","otherk_mygp_ctrl_prop","otherk_mygp_case_prop","my0_intactk_StartL_mean","fwd_intactk_count","rev_intactk_count")

#extract the unique kmer
myk4plot<-unique(myflk_behave_pheno$kmer)

for (j in 1:length(myk4plot)){ #open bracket for looping through each kmer
  #for (j in 1:10){ #open bracket for looping through each kmer
  
  mykmer<-myk4plot[j] 
  #mykmer<-"kmer964" 
  
  #print(mykmer)
  
  #select the rows referring to the kmer
  mytable<-myflk_behave_pheno[which(myflk_behave_pheno$kmer==mykmer),]
  
  #making the output matrix
  myout<-matrix(0,1,9)
  colnames(myout)<-c("kmer","event_sum","flk_behaviour","my0_intactk_sum","my1_intactk_sum","otherk","my0_otherk_sum","my1_otherk_sum","event")
  myout<-as.data.frame(myout)
  myout$kmer<-mykmer

#making the output matrix for short all out
  myshortout<-matrix(0,1,8)
  colnames(myshortout)<-c("kmer","intactk_mygp_ctrl_prop","intactk_mygp_case_prop","otherk_mygp_ctrl_prop","otherk_mygp_case_prop","my0_intactk_StartL_mean","fwd_intactk_count","rev_intactk_count")
  myshortout<-as.data.frame(myshortout)
  myshortout$kmer<-mykmer  
  
  #get the total number of cases and controls
  ctrl_count<-length(which(mytable$case_control=="0"))
  case_count<-length(which(mytable$case_control=="1"))
  
  #count the proportion of cases and controls in mybehave column
  mytable$mybehave<-as.character(mytable$mybehave) 
  mycat<-unique(mytable$mybehave)
  
  myout$event_sum<-paste(mycat,collapse=":") #fill in the table
  
  mysum_str<-""  #pasting different behaviours' count and proportion into one string 
  
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
    
    #get the coordinate summary statistics (based on StartL and StartR) of 
	
	if(any(c(mygp_ctrl_prop, mygp_case_prop)>0.2) & (mygp=="intact_k")){ #check that this behaviour is not the minority ones
	
	my0_intactk_StartL_stat<-paste(round(summary(mytable[which(mytable$mybehave==mygp & mytable$case_control=="0"),"StartL"]),0),collapse=" ")
	my0_intactk_StartL_SD<-round(sd(mytable[which(mytable$mybehave==mygp & mytable$case_control=="0"),"StartL"]),0)
	
	my0_intactk_StartR_stat<-paste(round(summary(mytable[which(mytable$mybehave==mygp & mytable$case_control=="0"),"StartR"]),0),collapse=" ")
	my0_intactk_StartR_SD<-round(sd(mytable[which(mytable$mybehave==mygp & mytable$case_control=="0"),"StartR"]),0)
	
	my0_intactk_flkdist<-paste(round(summary(mytable[which(mytable$mybehave==mygp & mytable$case_control=="0"),"flk_dist"]),0),collapse=" ")
	
	myout$my0_intactk_sum<-paste("StartL_stat:",my0_intactk_StartL_stat," | StartL_sd:",my0_intactk_StartL_SD, " | StartR_stat:", my0_intactk_StartR_stat, " | StartR_sd:", my0_intactk_StartR_SD," | flk_dist_stat:", my0_intactk_flkdist,collapse=" ")
	
	my1_intactk_StartL_stat<-paste(round(summary(mytable[which(mytable$mybehave==mygp & mytable$case_control=="1"),"StartL"]),0),collapse=" ")
	my1_intactk_StartL_SD<-round(sd(mytable[which(mytable$mybehave==mygp & mytable$case_control=="1"),"StartL"]),0)
	
	my1_intactk_StartR_stat<-paste(round(summary(mytable[which(mytable$mybehave==mygp & mytable$case_control=="1"),"StartR"]),0),collapse=" ")
	my1_intactk_StartR_SD<-round(sd(mytable[which(mytable$mybehave==mygp & mytable$case_control=="1"),"StartR"]),0)
	
	my1_intactk_flkdist<-paste(round(summary(mytable[which(mytable$mybehave==mygp & mytable$case_control=="1"),"flk_dist"]),0),collapse=" ")
	
	myout$my1_intactk_sum<-paste("StartL_stat:",my1_intactk_StartL_stat," | StartL_sd:",my1_intactk_StartL_SD, " | StartR_stat:", my1_intactk_StartR_stat, " | StartR_sd:", my1_intactk_StartR_SD," | flk_dist_stat:", my1_intactk_flkdist,collapse=" ")

        #for myshortout
	myshortout$intactk_mygp_ctrl_prop<-mygp_ctrl_prop
	myshortout$intactk_mygp_case_prop<-mygp_case_prop
        myshortout$my0_intactk_StartL_mean<-signif(as.numeric(mean(mytable[which(mytable$mybehave==mygp & mytable$case_control=="0"),"StartL"])),digits = as.numeric(as.character(opt$dedupk)))	
	}
	
	if(any(c(mygp_ctrl_prop, mygp_case_prop)>0.2) & (mygp!="intact_k")){ #check that this behaviour is not the minority ones
	
	myout$otherk<-as.character(mygp)
	
	my0_otherk_StartL_stat<-paste(round(summary(mytable[which(mytable$mybehave==mygp & mytable$case_control=="0"),"StartL"]),0),collapse=" ")
	my0_otherk_StartL_SD<-round(sd(mytable[which(mytable$mybehave==mygp & mytable$case_control=="0"),"StartL"]),0)
	
	my0_otherk_StartR_stat<-paste(round(summary(mytable[which(mytable$mybehave==mygp & mytable$case_control=="0"),"StartR"]),0),collapse=" ")
	my0_otherk_StartR_SD<-round(sd(mytable[which(mytable$mybehave==mygp & mytable$case_control=="0"),"StartR"]),0)
	
	my0_otherk_flkdist<-paste(round(summary(mytable[which(mytable$mybehave==mygp & mytable$case_control=="0"),"flk_dist"]),0),collapse=" ")
	
	myout$my0_otherk_sum<-paste("StartL_stat:",my0_otherk_StartL_stat," | StartL_sd:",my0_otherk_StartL_SD, " | StartR_stat:", my0_otherk_StartR_stat, " | StartR_sd:", my0_otherk_StartR_SD," | flk_dist_stat:", my0_otherk_flkdist,collapse=" ")
	
	my1_otherk_StartL_stat<-paste(round(summary(mytable[which(mytable$mybehave==mygp & mytable$case_control=="1"),"StartL"]),0),collapse=" ")
	my1_otherk_StartL_SD<-round(sd(mytable[which(mytable$mybehave==mygp & mytable$case_control=="1"),"StartL"]),0)
	
	my1_otherk_StartR_stat<-paste(round(summary(mytable[which(mytable$mybehave==mygp & mytable$case_control=="1"),"StartR"]),0),collapse=" ")
	my1_otherk_StartR_SD<-round(sd(mytable[which(mytable$mybehave==mygp & mytable$case_control=="1"),"StartR"]),0)
	
	my1_otherk_flkdist<-paste(round(summary(mytable[which(mytable$mybehave==mygp & mytable$case_control=="1"),"flk_dist"]),0),collapse=" ")
	
	myout$my1_otherk_sum<-paste("StartL_stat:",my1_otherk_StartL_stat," | StartL_sd:",my1_otherk_StartL_SD, " | StartR_stat:", my1_otherk_StartR_stat, " | StartR_sd:", my1_otherk_StartR_SD," | flk_dist_stat:", my1_otherk_flkdist,collapse=" ")
	
        #for myshortout
	myshortout$otherk_mygp_ctrl_prop<-mygp_ctrl_prop
	myshortout$otherk_mygp_case_prop<-mygp_case_prop
	}
    
    } #close bracket for looping through each behaviour
    
    myout$flk_behaviour<-mysum_str   #fill in the table with the behaviour summary
    
 if(myout$otherk%in%c("mv_away","swp_flk")){
myout$event<-"translocation"}
if(myout$otherk%in%c("mv&flp","swp&flp")){
myout$event<-"inversion"}

#finally count the number of few_intactk and rev_intactk for this kmer
myshortout$fwd_intactk_count<-length(which(mytable$intactk_orien=="intactk_fwd"))
myshortout$rev_intactk_count<-length(which(mytable$intactk_orien=="intactk_rev"))

  myall_out<-rbind(myall_out,myout)
  myshort_allout<-rbind(myshort_allout,myshortout)

} #close bracket for looping through each kmer

#keeping only the unique row for myshort_allout
myshort_allout$label<-apply( myshort_allout[,2:8] , 1 , paste , collapse = "-" )
myshort_allout_uniq<-myshort_allout[!duplicated(myshort_allout$label),]
myshort_allout_uniq<-myshort_allout_uniq[-1,]

write.table(myall_out,file=paste(opt$outdir,"mysplitk_out.txt",sep="/"),quote=F,row.names = F,col.names = T,sep="\t")
write.table(myshort_allout_uniq,file=paste(opt$outdir,"myshort_splitk_out_uniq.txt",sep="/"),quote=F,row.names = F,col.names = T,sep="\t")
}else{
print("no split k")
}

print("make summary for intactk")

#now process table of kmers with flank behaviour of "intact_k" only, also remove "StartL" and "EndR" column to have less to process
#the theory is these intact kmers is flagged as significantly associated in GWAS because it is within the "inverted" genomic region
#aim to find association between fwd_k and rev_k and "0" and "1", i.e. 
#the proportion of "0" genomes with fwd_k
#the proportion of "0" genomes with rev_k
#the proportion of "1" genomes with fwd_k
#the proportion of "1" genomes with rev_k
#let's say fwd_k is associated with "0", then if the SD of the genomic position is 
#very small (i.e. suggesting the same position), then plot the medium genomic 
#position of where the fwd_k is mapped in "0" genomes

#specify the number of genomes
numgen<-length(mygen)
#numgen<-47

myallk<-unique(myall_behave$kmer)

#list of kmers with flank behaviour == "intact_k"
myk_intactbebave<-myallk[which(myallk%!in%myk_otherbebave)]

if(length(myk_intactbebave)>0){
myintactk_only_tab<-myall_behave[which(myall_behave$kmer%in%myk_intactbebave),-c(4,7)]

#add the new columns for storing kmer orientation information
myintactk_only_tab$k_orien<-NA

#define the forward intact k
myfwdk_k_coor<-which(myintactk_only_tab$EndL<myintactk_only_tab$StartR)
myintactk_only_tab[myfwdk_k_coor,"k_orien"]<-"fwd_k"

myrevk_k_coor<-which(myintactk_only_tab$EndL>myintactk_only_tab$StartR)
myintactk_only_tab[myrevk_k_coor,"k_orien"]<-"rev_k"

#then check for each kmer if fwd_k is found in most "0" pheno and most rev_k is found in most "1" pheno, and vice versa
#also check if certain location is associated with "0"/"1" pheno

#first check if the rows for each kmer (one row per genome) is the same as the total number of genomes used in GWAS
myfreq<-as.data.frame(table(myintactk_only_tab$kmer))
table(myfreq$Freq)
#    0    47 
#  504 10220 

#keep those kmers with blast hit in all genomes only (list of kmers)
myk4paint<-myfreq[which(myfreq$Freq==numgen),"Var1"]   #change genome number here

#making intact_k summary output table for all intact kmer
myintactk_out<-matrix(0,1,21)

colnames(myintactk_out)<-c("kmer","kmer_behaviour","flk_dist","fwdk_gen_count","revk_gen_count","fwdk_0gen_prop","revk_0gen_prop","fwdk_1gen_prop","revk_1gen_prop","fwdk_0gen_count","revk_0gen_count","fwdk_1gen_count","revk_1gen_count","fwdk_0gen_med","fwdk_0gen_sd","revk_0gen_med","revk_0gen_sd","fwdk_1gen_med","fwdk_1gen_sd","revk_1gen_med","revk_1gen_sd")

for (i in 1:length(myk4paint)){

mykmer<-myk4paint[i]
#mykmer<-"kmer9997"

#print(as.character(mykmer))

#get the number of "0" genomes and "1" genomes by looking at the first kmer
myzero<-length(which(myintactk_only_tab$kmer==myk4paint[1] & myintactk_only_tab$case_control==0))
myone<-length(which(myintactk_only_tab$kmer==myk4paint[1] & myintactk_only_tab$case_control==1))

#get the rows of the kmer
mysub<-myintactk_only_tab[which(myintactk_only_tab$kmer==mykmer),]

#my unique behave 
mybehave<-toString(unique(mysub$mybehave))

#my unique flk dist 
myflk_dist<-toString(unique(mysub$flk_dist))

#find out the number of genome with fwd_k
myfwd_k_count<-length((which(mysub$k_orien=="fwd_k")))

#find out the number of genome with rev_k
myrev_k_count<-length((which(mysub$k_orien=="rev_k")))

#find out if fwd_k and rev_k is associated with "0" and "1" genomes
#for this kmer, the proportion of "0" genomes with fwd_k
my0_fwdk_prop<-round(nrow(mysub[which(mysub$case_control==0 & mysub$k_orien=="fwd_k"),])/myzero,3)
#for this kmer, the proportion of "0" genomes with rev_k
my0_revk_prop<-round(nrow(mysub[which(mysub$case_control==0 & mysub$k_orien=="rev_k"),])/myzero,3)
#for this kmer, the proportion of "1" genomes with fwd_k
my1_fwdk_prop<-round(nrow(mysub[which(mysub$case_control==1 & mysub$k_orien=="fwd_k"),])/myone,3)
#for this kmer, the proportion of "1" genomes with rev_k
my1_revk_prop<-round(nrow(mysub[which(mysub$case_control==1 & mysub$k_orien=="rev_k"),])/myone,3)

#the number of "0" genomes with fwd_k
my0_fwdk_count<-nrow(mysub[which(mysub$case_control==0 & mysub$k_orien=="fwd_k"),])
#fthe number of "0" genomes with rev_k
my0_revk_count<-nrow(mysub[which(mysub$case_control==0 & mysub$k_orien=="rev_k"),])
#the number of "1" genomes with fwd_k
my1_fwdk_count<-nrow(mysub[which(mysub$case_control==1 & mysub$k_orien=="fwd_k"),])
#the number of "1" genomes with rev_k
my1_revk_count<-nrow(mysub[which(mysub$case_control==1 & mysub$k_orien=="rev_k"),])

#describe the median and SD genomic position of fwd_k + "0" genomes 
my0_fwdk_medium<-median(mysub[which(mysub$case_control==0 & mysub$k_orien=="fwd_k"),"EndL"])
my0_fwdk_sd<-sd(mysub[which(mysub$case_control==0 & mysub$k_orien=="fwd_k"),"EndL"])

#describe the median and SD genomic position of rev_k + "0" genomes 
my0_revk_medium<-median(mysub[which(mysub$case_control==0 & mysub$k_orien=="rev_k"),"EndL"])
my0_revk_sd<-sd(mysub[which(mysub$case_control==0 & mysub$k_orien=="rev_k"),"EndL"])

#describe the median and SD genomic position of fwd_k + "1" genomes 
my1_fwdk_medium<-median(mysub[which(mysub$case_control==1 & mysub$k_orien=="fwd_k"),"EndL"])
my1_fwdk_sd<-sd(mysub[which(mysub$case_control==1 & mysub$k_orien=="fwd_k"),"EndL"])

#describe the median and SD genomic position of rev_k + "1" genomes 
my1_revk_medium<-median(mysub[which(mysub$case_control==1 & mysub$k_orien=="rev_k"),"EndL"])
my1_revk_sd<-sd(mysub[which(mysub$case_control==1 & mysub$k_orien=="rev_k"),"EndL"])

myrowout<-c(as.character(mykmer),mybehave,myflk_dist,myfwd_k_count,myrev_k_count,my0_fwdk_prop,my0_revk_prop,my1_fwdk_prop,my1_revk_prop,my0_fwdk_count,my0_revk_count,my1_fwdk_count,my1_revk_count,my0_fwdk_medium,my0_fwdk_sd,my0_revk_medium,my0_revk_sd,my1_fwdk_medium,my1_fwdk_sd,my1_revk_medium,my1_revk_sd)

myintactk_out<-rbind(myintactk_out,myrowout)

}

write.table(myintactk_out,file=paste(opt$outdir,"myintactkwithN_out.txt",sep="/"),quote=F,row.names = F,col.names = T,sep="\t")
}else{
print("no intact k")
}

}  #close bracket for process good kmers
