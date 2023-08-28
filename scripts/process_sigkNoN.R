#usage: Rscript process_sigkNoN.R --intactk  --gen --pheno  --outdir
#Rscript process_sigkNoN.R --pheno prn_status_pheno.txt \
#--outdir /home/ubuntu/Dorothy/genome_rearrangement/PRN_468_tree_merge200GWAS_noNsigk/test_dir_28Mar23

library("optparse")

option_list = list(
  make_option("--pheno", type="character", default=NULL,
              help="kmer to plot", metavar="character"),
  make_option("--outdir", type="character", default=NULL,
              help="kmer to plot", metavar="character")
) 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); 

#headers from the blast output
#query, subject, identity, alig_len, mismatches, gap, qstart, qend, sStart, sEnd, evalue, bitscore

###blast output quality check
#all kmer that should show at least one blast hit

#setwd(opt$outdir)

#load in phenotype file
myphenofile<-read.table(opt$pheno,header=F)
#myphenofile<-read.table("prn_status_pheno.txt",header=F)  #check that there should be no header in phenotype file

#load in the length of the sig kmers
myflk_coor<-read.delim(paste(opt$outdir,"kmernoN_length.txt",sep="/"),header=F,sep="_")
colnames(myflk_coor)<-c("kmer","kmer_len")
head(myflk_coor)

#load in the blast output file
mytable<-read.table(paste(opt$outdir,"mynoN_out.txt",sep="/"), header=F)
colnames(mytable)<-c("query","subject","identity","alig_len","mismatches","gap","qstart","qend","sStart","sEnd","evalue","bitscore")

mykmer<-as.character(unique(mytable$query))  #get the list of kmers with blast output
mygen<-as.character(unique(myphenofile$V1))  #get the list of the genomes from the pheno file

#set the output for the rows of kmers with absence of >5% genomes
abs_gen_k<-c()

#set the output for the rows with multiple blast hits in flank
multi_hit_k<-c()

#set the output for the rows with alignment length issue
align_len_k<-c()

#set the output for the rows with identity and  e value issue
ID_E_issue_k<-c()

#looping through kmers
for (i in 1:length(mykmer)){

myalignlen<-myflk_coor$kmer_len[which(as.character(myflk_coor$kmer)==as.character(mykmer[i]))]
#print(myalignlen)

#print(mykmer[i])
mykrow<-mytable[which(mytable$query==as.character(mykmer[i])),]
#mykrow<-mytable[which(mytable$query=="kmer51"),]

#blast hit must contain >=95% of the genomes used
if(length(unique(mykrow$subject))<(0.95*length(mygen))){
abs_gen_k<-c(abs_gen_k,mykmer[i])
}

#each genome must appear once only
myfreqtable<-data.frame(table(mykrow$subject))
if(any(myfreqtable$Freq>1)){
multi_hit_k<-c(multi_hit_k,mykmer[i])
}

#each blast match should be at least 90% in length of the flank length
mykrow$alignlen<-NA
mykrow$alignlen<-mykrow$alig_len/myalignlen
#print(mykrow$alignlen)
  if(any(mykrow$alignlen<0.9)){
    align_len_k<-c(align_len_k,mykmer[i])
  }

#use percentage identity and E value filter
if(any(mykrow$identity<95) | any(mykrow$evalue>10e-10)){
  ID_E_issue_k<-c(ID_E_issue_k,mykmer[i])
  }

}#close the for loop

#make summary for kmer quality control
myfilterout<-matrix(0,length(mykmer),4)
rownames(myfilterout)<-mykmer
colnames(myfilterout)<-c("abs_gen_k","multi_hit_k","align_len_k","ID_E_issue_k")

myfilterout<-as.data.frame(myfilterout)
#myfilterout

#output the rows of kmers with quality issues

if (length(unique(abs_gen_k))>0){
print(paste("abs_gen_k has ",length(unique(abs_gen_k))," kmers",sep=""))
myfilterout[abs_gen_k,"abs_gen_k"]<-"yes"
my_abs_gen_k<-mytable[which(mytable$query%in%unique(abs_gen_k)),]
write.table(my_abs_gen_k,file=paste(opt$outdir,"kmer_with_missinggenomes_noN.txt",sep="/"),quote=F,row.names=F,col.names = T,sep="\t")
}else{
print("no abs_gen_k")
}


if (length(unique(multi_hit_k))>0){
print(paste("multi_hit_k has ",length(unique(multi_hit_k))," kmers",sep=""))
myfilterout[multi_hit_k,"multi_hit_k"]<-"yes"
my_multi_hit_k<-mytable[which(mytable$query%in%unique(multi_hit_k)),]
write.table(my_multi_hit_k,file=paste(opt$outdir,"kmer_with_multi_hits_noN.txt",sep="/"),quote=F,row.names=F,col.names = T,sep="\t")
}else{
print("no multi_hit_k")
}

if (length(unique(align_len_k))>0){
print(paste("align_len_k has ",length(unique(align_len_k))," kmers",sep=""))
myfilterout[align_len_k,"align_len_k"]<-"yes"
my_align_len_k<-mytable[which(mytable$query%in%unique(align_len_k)),]
write.table(my_align_len_k,file=paste(opt$outdir,"kmer_with_align_len_noN.txt",sep="/"),quote=F,row.names=F,col.names = T,sep="\t")
}else{
print("no align_len_k")
}

if (length(unique(ID_E_issue_k))>0){
print(paste("ID_E_issue_k has ",length(unique(ID_E_issue_k))," kmers",sep=""))
myfilterout[ID_E_issue_k,"ID_E_issue_k"]<-"yes"
my_ID_E_issue_k<-mytable[which(mytable$query%in%unique(ID_E_issue_k)),]
write.table(my_ID_E_issue_k,file=paste(opt$outdir,"kmer_with_ID_E_issue_noN.txt",sep="/"),quote=F,row.names=F,col.names = T,sep="\t")
}else{
print("no ID_E_issue_k")
}

write.table(myfilterout,file=paste(opt$outdir,"filterk_out_summary_noN.txt",sep="/"),quote=F,row.names = T,col.names =T,sep="\t")


#output the good kmers for further processing
mybadk<-unique(c(abs_gen_k,multi_hit_k,align_len_k,ID_E_issue_k))

mygoodk<-mykmer[which(!is.element(mykmer,mybadk))]

myprocess<-mytable[which(mytable$query%in%mygoodk),]

#the myprocess table should refere to kmers that are present in all genomes with both flanks; the flanks are also fully aligned with no SNPs nor gaps, and the flanks show unique blast hit in each genome 
write.table(myprocess,file=paste(opt$outdir,"rows_for_process_NoN.txt",sep="/"),quote=F,row.names = F,col.names = T,sep="\t")

#read in myprocess table
myprocess<-read.table(paste(opt$outdir,"rows_for_process_NoN.txt",sep="/"),sep="\t",header=T)

#determining orientation for each blast hit
myprocess$k_orien<-0
myprocess[which(myprocess$sEnd>myprocess$sStart),"k_orien"]<-"fwd_k"
myprocess[which(myprocess$sStart>myprocess$sEnd),"k_orien"]<-"rev_k"

mystartendLR_out<-myprocess[,c("query","subject","sStart","sEnd","k_orien")]
colnames(mystartendLR_out)<-c("kmer","genome","sStart","sEnd","k_orien")

#merge in the phenotype information
myflk_behave_pheno<-merge(mystartendLR_out,myphenofile,by.x="genome",by.y="V1",)

myflk_behave_pheno<-myflk_behave_pheno[order(myflk_behave_pheno$kmer, decreasing=T),]
colnames(myflk_behave_pheno)[6]<-"case_control"

write.table(myflk_behave_pheno,file=paste(opt$outdir,"myflk_behave_pheno_NoN.txt",sep="/"),quote=F,row.names = F,col.names = T,sep="\t")

#read in myflk_behave_pheno.txt
myintactk_only_tab<-read.table(paste(opt$outdir,"myflk_behave_pheno_NoN.txt",sep="/"),sep="\t",header=T)
colnames(myintactk_only_tab)[5]<-"k_orien"

#now process table of kmers with flank behaviour of "intact_k" only, 
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
#numgen<-468

#then check for each kmer if fwd_k is found in most "0" pheno and most rev_k is found in most "1" pheno, and vice versa
#also check if certain location is associated with "0"/"1" pheno

#first check if the rows for each kmer (one row per genome) is the same as the total number of genomes used in GWAS
myfreq<-as.data.frame(table(myintactk_only_tab$kmer))
table(myfreq$Freq)
#    0    47 
#  504 10220 

#keep those kmers with blast hit in all genomes only (list of kmers), optional
#myk4paint<-myfreq[which(myfreq$Freq==numgen),"Var1"]   #change genome number here

#get the list of intact kmers for process
myk4paint<-myfreq[,"Var1"] 

#making intact_k summary output table for all intact kmer
myintactk_out<-matrix(0,1,21)

colnames(myintactk_out)<-c("kmer","kmer_behaviour","flk_dist","fwdk_gen_count","revk_gen_count","fwdk_0gen_prop","revk_0gen_prop","fwdk_1gen_prop","revk_1gen_prop","fwdk_0gen_count","revk_0gen_count","fwdk_1gen_count","revk_1gen_count","fwdk_0gen_med","fwdk_0gen_sd","revk_0gen_med","revk_0gen_sd","fwdk_1gen_med","fwdk_1gen_sd","revk_1gen_med","revk_1gen_sd")

for (i in 1:length(myk4paint)){

mykmer<-myk4paint[i]
#mykmer<-"kmer9997"

#print(as.character(mykmer))

mybehave<-"intactk_noN"

#get the number of "0" genomes and "1" genomes by looking at the first kmer
myzero<-length(which(myintactk_only_tab$kmer==myk4paint[1] & myintactk_only_tab$case_control==0))
myone<-length(which(myintactk_only_tab$kmer==myk4paint[1] & myintactk_only_tab$case_control==1))

#get the rows of the kmer
mysub<-myintactk_only_tab[which(myintactk_only_tab$kmer==mykmer),]

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
my0_fwdk_medium<-median(mysub[which(mysub$case_control==0 & mysub$k_orien=="fwd_k"),"sStart"])
my0_fwdk_sd<-sd(mysub[which(mysub$case_control==0 & mysub$k_orien=="fwd_k"),"sStart"])

#describe the median and SD genomic position of rev_k + "0" genomes 
my0_revk_medium<-median(mysub[which(mysub$case_control==0 & mysub$k_orien=="rev_k"),"sStart"])
my0_revk_sd<-sd(mysub[which(mysub$case_control==0 & mysub$k_orien=="rev_k"),"sStart"])

#describe the median and SD genomic position of fwd_k + "1" genomes 
my1_fwdk_medium<-median(mysub[which(mysub$case_control==1 & mysub$k_orien=="fwd_k"),"sStart"])
my1_fwdk_sd<-sd(mysub[which(mysub$case_control==1 & mysub$k_orien=="fwd_k"),"sStart"])

#describe the median and SD genomic position of rev_k + "1" genomes 
my1_revk_medium<-median(mysub[which(mysub$case_control==1 & mysub$k_orien=="rev_k"),"sStart"])
my1_revk_sd<-sd(mysub[which(mysub$case_control==1 & mysub$k_orien=="rev_k"),"sStart"])

myrowout<-c(as.character(mykmer),mybehave,"NA",myfwd_k_count,myrev_k_count,my0_fwdk_prop,my0_revk_prop,my1_fwdk_prop,my1_revk_prop,my0_fwdk_count,my0_revk_count,my1_fwdk_count,my1_revk_count,my0_fwdk_medium,my0_fwdk_sd,my0_revk_medium,my0_revk_sd,my1_fwdk_medium,my1_fwdk_sd,my1_revk_medium,my1_revk_sd)

myintactk_out<-rbind(myintactk_out,myrowout)

}

myintactk_out<-myintactk_out[-1,]

write.table(myintactk_out,file=paste(opt$outdir,"myNoNintactk_out.txt",sep="/"),quote=F,row.names = F,col.names = T,sep="\t")
