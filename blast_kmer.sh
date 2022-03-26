###main inputs for the pipeline###
#total number of genomes (input of the script)
gen_num<-49

#kmers'length  (input of the script)
k_len=200

####################################################

#each flank needs to be at least 30bp long for the kmer to be processed
#make phenotype file for pipeline
sed '1d' phenotypes.tsv > case_control.txt

#input should be a kmer multifasta file

#extract the sig kmer lines from the output file
#sig. kmer without bad-chisq with N
awk '{ if ($4 <= 4.13E-04) { print } }' chromstruc_clus1clus2_GWAS_nopopctrl | grep "N" > allsigkmer_withN.txt

#extract the sig kmer sequences

awk 'NR!=1{print $1}' /home/ubuntu/Dorothy/USAgenomes_GWAS/chromstruc_clus1clus2_GWAS/nopopctrl/maf0.05_replremainIS/allsigkmer_withN.txt > sig_kmer_list.txt

#making the headerfile
number=$(cat /home/ubuntu/Dorothy/USAgenomes_GWAS/chromstruc_clus1clus2_GWAS/nopopctrl/maf0.05_replremainIS/allsigkmer_withN.txt | wc -l)

START=1
let "END=$number-1"  #not including the header line
 
for (( c=$START; c<=$END; c++ ))
do
	echo ">kmer""$c " >> header.txt
done

#create a multifasta file of the associated kmers by joining the two files alternatively
paste -d \\n header.txt sig_kmer_list.txt > allsig_kmer.fasta

##start of the main loop##

#set output directory
out_path=/home/ubuntu/Dorothy/USAgenomes_GWAS/chromstruc_clus1clus2_GWAS/nopopctrl/maf0.05_replremainIS  #set the output directory of the blast output

#get the flank start and end coordinates of the sig kmers, output file name: flank_coor.txt
python3 /home/ubuntu/Dorothy/USAgenomes_GWAS/scripts_pipeline/extract_flank_coor.py --input allsig_kmer.fasta

#remove kmers with short flanks and extract kmers for blasting, by looping through the header file
while read mykmer
do
mykmer=$(echo ${mykmer} | sed 's/>//')  #remove ">" from the header line
#mykmer="kmer1" # specify which kmer being process
echo ${mykmer}

#filter the kmers that contain each flanks of 30bp
leftflankend=$(grep ${mykmer}_ flank_coor.txt | cut -d "_" -f2)
rightflankstart=$(grep ${mykmer}_ flank_coor.txt | cut -d "_" -f3)
kmerlen=$(grep ${mykmer}_ flank_coor.txt | cut -d "_" -f4)
rightlen="$((${kmerlen}-${rightflankstart}+1))"       

if [ ${leftflankend} -lt 30 ] || [ ${rightlen} -lt 30 ]; then   #filter out the kmer with at least one flank of <30bp
grep ${mykmer}_ flank_coor.txt >> $out_path/kmer_flanktooshort_flkcoor.txt  #store the flank coordinates of kmer with too short flank
echo ${mykmer} >> $out_path/kmer_flanktooshort_4rm.txt   #storing the kmer with flank being too short for remove
fi

done  < /home/ubuntu/Dorothy/USAgenomes_GWAS/chromstruc_clus1clus2_GWAS/nopopctrl/maf0.05_replremainIS/header.txt #closing for looping through sig kmers header lines]


awk 'BEGIN{while((getline<"kmer_flanktooshort_4rm.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' allsig_kmer_withN.fasta > kmer_forblast.fasta #remove the kmer with flanks being too short from the muktifasta file for blasting


#make a multifasta file of the genomes for blasting
#moving the relevant genomes into one directory
mkdir /home/ubuntu/Dorothy/B.pertussis_573genomes_NCBI/Weigand_USA_genome/107usa_samplelist
while read line
do
mv /home/ubuntu/Dorothy/B.pertussis_573genomes_NCBI/Weigand_USA_genome/${line}.fna.gz /home/ubuntu/Dorothy/B.pertussis_573genomes_NCBI/Weigand_USA_genome/107usa_samplelist
done < /home/ubuntu/Dorothy/USAgenomes_GWAS/107_typeGWAS/107usa_samplelist.txt

#renaming the genome seqID as the isolate name
#(python)
import argparse
from Bio import SeqIO

filepath="/home/ubuntu/Dorothy/USAgenomes_GWAS/107_typeGWAS/107usa_samplelist.txt"
with open(filepath) as fp: 
      count=0
      sample = open(filepath,'r').read().split('\n') 
      del sample[-1]
      print(sample)

for i in range(0,len(sample)):
	print(str(sample[i]))
	fasta_sequences = SeqIO.parse(str(sample[i])+".fna",'fasta') #this input is the output of the above script
	for record in fasta_sequences:
		record.description=str(record.id)+" "+str(record.description)
		record.id=str(sample[i])
		SeqIO.write(record, str(sample[i])+"_rename.fna", "fasta")


#blasting, subjects are genomes in multifasta file, genome names should be in the headers 
blastn -query kmer_forblast.fasta -subject /home/ubuntu/Dorothy/B.pertussis_573genomes_NCBI/Weigand_USA_genome/49_clusterGWAS_genlist.fasta -outfmt 6 -out myout.txt
  
  
####################################################

#headers from the blast output
#query, subject, identity, alig_len, mismatches, gap, qstart, qend, sStart, sEnd, evalue, bitscore

###kmer hits quality checks###

#(R)

#make sure the input file kmer_blastcoor.txt is in the working directory
mytable<-read.table("myout.txt", header=F)
colnames(mytable)<-c("query","subject","identity","alig_len","mismatches","gap","qstart","qend","sStart","sEnd","evalue","bitscore")

#load in the phenotype file 
myphenofile<-read.delim("/home/ubuntu/Dorothy/USAgenomes_GWAS/chromstruc_clus1clus2_GWAS/case_control.txt", header=F)

#load in the flank start and end coordinates of the sig kmers
myflk_coor<-read.delim("flank_coor.txt",header=F,sep="_")
colnames(myflk_coor)<-c("kmer","leftflankend","rightflankstart","kmer_len")

mykmer<-as.character(unique(mytable$query))  #get the list of kmers with blast output, all kmer that should show at least one blast hit
mygen<-as.character(unique(myphenofile$V1))  #get the list of the genomes from the pheno file

#set the output for the rows of kmers with deletions
del_k<-c()

#set the output for the rows with multiple blast hits in flank 
multi_hit_k<-c()

#set the output for the rows with incomplete flank alignment
alignlen_issue_k<-c()

#set the output for the rows with SNPs and gaps
SNPgap_k<-c()

#create matrix for storing the rows referring to samples with presence of both flanks
myprocess<-matrix(0,0,ncol(mytable))

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
if (length(unique(myk_blastcoor))==4 & (any(is.element(myk_blastcoor,myflk_coor_k)!=T | any(is.element(myflk_coor_k,myk_blastcoor)!=T)))){  #output the kmers with incomplete flank alignment
 alignlen_issue_k<-c(alignlen_issue_k,mykmer[i])
 count=count+1
 } 
 if (any(mysnpgap!=0)) {  #output the kmers with SNP or gaps in flanks
 SNPgap_k<-c(SNPgap_k,mykmer[i])
 count=count+1
 } 
 
 } #closing for loop for genomes
 
 if (count==0){
 myprocess<-rbind(myprocess,myk_row)
 }
 
 } #closing for loop for kmers
 

#the myprocess table should refere to kmers that are present in all genomes with both flanks; the flanks are also fully aligned with no SNPs nor gaps, and the flanks show unique blast hit in each genome

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

write.table(myprocess,file="rows_for_process.txt",quote=F,row.names = F,col.names = T,sep="\t")



#sort by the samples
#sort -k2 -h kmer2_blastcoor.txt > kmer2_blastcoor_sort.txt  
