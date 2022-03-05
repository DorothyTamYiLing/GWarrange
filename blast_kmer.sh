#I am going to blast multifasta kmer and multifasta genomes files
#each flank needs to be at least 30bp long for the kmer to be processed

#creating the genome path file, IF the file is present, dont run it again!!
while read samples   #reading in each genome  
do 
echo /home/ubuntu/Dorothy/B.pertussis_573genomes_NCBI/Weigand_USA_genome/${samples}.fna >> /home/ubuntu/Dorothy/USAgenomes_GWAS/107_typeGWAS/all_path.txt
done < /home/ubuntu/Dorothy/USAgenomes_GWAS/107_typeGWAS/107usa_samplelist.txt

#input should be a kmer multifasta file

#extract the sig kmer lines from the output file
#sig. kmer without bad-chisq with N
awk '{ if ($4 <= 3.42E-04) { print } }' 107USA_typeGWAS_ISrepl_maf0.05_fix_nobad-chisq | grep "N" > allsigkmer_withN.txt

#extract the sig kmer sequences

awk 'NR!=1{print $1}' allsigkmer_withN.txt > sig_kmer_list.txt

#making the headerfile
number=$(cat allsigkmer_withN.txt | wc -l)

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
out_path=/home/ubuntu/Dorothy/USAgenomes_GWAS/chromstruc_clus1clus2_GWAS/nopopctrl/maf0.2/gen_rearr_pipeline  #set the output directory of the blast output

#extract each kmer sequence for processing by looping through the header file
while read mykmer
do
#mykmer=$(echo ${mykmer} | sed 's/>//')
#echo ${mykmer}
mykmer="kmer2" # specify which kmer being process
#done  < /home/ubuntu/Dorothy/USAgenomes_GWAS/chromstruc_clus1clus2_GWAS/nopopctrl/maf0.2/header.txt #closing for looping through sig kmers]

python3 /home/ubuntu/Dorothy/USAgenomes_GWAS/chromstruc_clus1clus2_GWAS/nopopctrl/maf0.2/scripts_pipeline/extract_kmer.py --allk /home/ubuntu/Dorothy/USAgenomes_GWAS/chromstruc_clus1clus2_GWAS/nopopctrl/maf0.2/allsig_kmer.fasta --kmer ${mykmer}

#creating the flank sequence files
python3 /home/ubuntu/Dorothy/USAgenomes_GWAS/chromstruc_clus1clus2_GWAS/nopopctrl/maf0.2/scripts_pipeline/slice_left_right_flk.py --input ${mykmer}.fasta

#looping from left flank to right flank [main loop starts here]
echo "mykmer genome genome_len flank flank_len genome_start genome_end flkstart flkend leftS/rightS leftE/rightE orientation identity evalue align_cov num_nt_aligned gaps num_SNP flk_presence" >> $out_path/${mykmer}_blastcoor.txt   

for x in {1..2}
do
if [ $x -eq 1 ] ; then
	myflankfile=${mykmer}_leftflank.fasta
	myflank="leftflank"
	echo processing left flank
fi
if [ $x -eq 2 ] ; then
	myflankfile=${mykmer}_rightflank.fasta
	myflank="rightflank"
	echo processing right flank
fi

#done

echo $myflankfile
echo $myflank


while read line   #reading in each genome  
do 
#echo $line    #line is the path including the genome file name

#line=/home/ubuntu/Dorothy/B.pertussis_573genomes_NCBI/Weigand_USA_genome/H698.fna.gz   #set the genome to be blasted
    gen_file=$(basename "${line}")  #get the genome fasta file name only
    #echo ${gen_file}
    gen_nm=${gen_file/.fna}  #get the genome name
    echo ${gen_nm}
    
    #unzip the genome left and right flank blast
    #gunzip ${line}  #unzip the genome fasta file
    
    #blast left flank
    out_file=${gen_nm}_${mykmer/.fasta}_${myflank}_output.txt  #set the output file name
    echo $out_file
    blastn -query ${line} -subject ${myflankfile} -out $out_path/$out_file
    
    #gzip the genome
    #gzip ${line/.gz}
    
    #detecting no hit found case
    if grep "No hits found" $out_path/$out_file; then
       echo flank shows no hit in the genome
       echo ${mykmer/.fasta} ${gen_nm} genome_len ${myflank} NA NA NA NA NA NA NA NA NA NA NA NA NA NA no_hit >> $out_path/${mykmer}_blastcoor.txt
    else
       csplit -f file -n 1 -k $out_path/$out_file /"bits"/ {*} #split the file with the occurrences of the string "bits", the output files will be replaced in the loop  
       file_count=$(ls file* | wc -l ) 
       let "number_hits=$file_count - 1"  #work
       if [ $file_count -gt 2 ]; then  #the flank should only show one blast hit, otherwise save into another file for investigation
           echo flank shows multiple hits in the genome
           #add lines here to output the multiple hit information (a for loop), start with for i in $(seq 1 $number_hits) 
           echo ${mykmer/.fasta} ${myflank} >> $out_path/flank_multiple_hit.txt
           echo ${mykmer/.fasta} ${gen_nm} genome_len ${myflank} NA NA NA NA NA NA NA NA NA NA NA NA NA NA multi_hit >> $out_path/${mykmer}_blastcoor.txt
       else  #if flank show one hit only
           #extract all the information from the output
           echo one hit in genome only
           gen_len=$(grep "Length=" file0 | head -1 | cut -d "=" -f2)  #length of the genome
           flk_len=$(grep "Length=" file0 | head -2 | tail -1 | cut -d "=" -f2)   #length of the flank
           identity=$(grep "Identities" file1 | head -n 1 | cut -d "(" -f2| cut -d "%" -f1)   #work
           flkstart=$(grep "Sbjct" file1 | head -n 1 | cut -d " " -f3)   #flank start
           flkend=$(grep "Sbjct" file1 | tail -n 1 | rev | cut -d " " -f1 | rev)  #flank end
           genstart=$(grep "Query" file1 | head -n 1 | cut -d " " -f3)   #genome start
           genend=$(grep "Query" file1 | tail -n 1 | rev | cut -d " " -f1 | rev)   #genome end  #work, grep the last field splitted by " "
           evalue=$(grep "Score" file1 | cut -d "=" -f3)
           alignment_coverage=$(grep "Identities" file1 | head -n 1 | cut -d "/" -f2 | cut -d "(" -f1) #work
           number_of_nt_aligned=$(grep "Identities" file1 | head -n 1 | cut -d "/" -f1 | cut -d "=" -f2) #work
           num_SNP=$(expr ${alignment_coverage}-${number_of_nt_aligned})
           gaps=$(grep "Identities" file1 | head -n 1 | cut -d "=" -f3 | cut -d "/" -f1) #work
           
           #defining if the flank is present
           if [ ${alignment_coverage} -eq ${flk_len} ] && [ ${num_SNP} -lt 1 ] && [ ${gaps} -lt 1 ]; then   #if the alignment coverage == flank length, and number of SNP or gap <1, define the flank to be present
               echo this flank is defined as present 
               presence="flank_present"
           else
               presence="flank_needs_investigation"
           fi 
           
           #defining the LeftS/RightS and LeftE/RightE, and check if the flank is flipped   , make sure genstart refer to the first position in flank
           if [ ${flkstart} -gt ${flkend} ]; then  #it means ${flkend}==1
                myS=${genend}   #therefore ${genend} should be the start of left flank
                myE=${genstart}
                orientation="reverse"
                echo ${myflank} is reversed
           else
                myS=${genstart}
                myE=${genend}
                orientation="forward"
                echo ${myflank} is forward
           fi 
           echo ${mykmer} ${gen_nm} ${gen_len} ${myflank} ${flk_len} ${genstart} ${genend} ${flkstart} ${flkend} ${myS} ${myE} ${orientation} ${identity} ${evalue} ${alignment_coverage} ${number_of_nt_aligned} ${gaps} ${num_SNP} ${presence} >> $out_path/${mykmer}_blastcoor.txt 
       fi
    fi 
    rm $out_path/$out_file
    rm file*
    
done < /home/ubuntu/Dorothy/USAgenomes_GWAS/107_typeGWAS/all_path.txt  #the genomes in this file are gzipped *gz

done  #closing for looping between left and right flank [main loop ends here]

#make phenotype file for pipeline
sed '1d' phenotypes.tsv > case_control.txt

Rscript /home/ubuntu/Dorothy/USAgenomes_GWAS/chromstruc_clus1clus2_GWAS/nopopctrl/maf0.2/scripts_pipeline/flank_behaviour.R 

rm *.fasta
rm kmer_behaviour_summary.txt
rm kmer_blastcoor.txt 

mv *_oneline_summary.txt all_oneline

done  < /home/ubuntu/Dorothy/USAgenomes_GWAS/chromstruc_clus1clus2_GWAS/nopopctrl/maf0.2/mytest.txt #closing for looping through sig kmers]

#sort by the samples
#sort -k2 -h kmer2_blastcoor.txt > kmer2_blastcoor_sort.txt  



