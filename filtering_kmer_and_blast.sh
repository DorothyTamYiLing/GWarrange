###filtering sig. kmers for blasting by keeping only the kmers that contain flanking sequences (both side) of at least $flnk_len bp in size

#$1 is k_input="/home/ubuntu/Dorothy/USAgenomes_GWAS/111_yearGWAS/k200/fix_model/allsig_kmer_withN.fasta"
#$2 is gen_input="/home/ubuntu/Dorothy/B.pertussis_573genomes_NCBI/Weigand_USA_genome/111_yearGWAS_genlist.fasta.gz"
#$3 is flnk_len=30

#get the flank start and end coordinates of the sig kmers, output file name: flank_coor.txt
python3 extract_flank_coor.py --input $1

#making the headerfile based on the number of significant kmers in k_input.fasta
k_num=$(grep ">" $1 | wc -l)

START=1
let "END=$k_num"  #not including the header line
 
for (( c=$START; c<=$END; c++ ))
do
	echo ">kmer""$c " >> header.txt
done

#remove kmers with short flanks and extract kmers for blasting, by looping through the header file
while read mykmer
do
mykmer=$(echo ${mykmer} | sed 's/>//')  #remove ">" from the header line
#mykmer="kmer1" # specify which kmer being process
#echo ${mykmer}

leftflankend=$(grep ${mykmer}_ flank_coor.txt | cut -d "_" -f2)
rightflankstart=$(grep ${mykmer}_ flank_coor.txt | cut -d "_" -f3)
kmerlen=$(grep ${mykmer}_ flank_coor.txt | cut -d "_" -f4)
rightlen="$((${kmerlen}-${rightflankstart}+1))"       

if [ ${leftflankend} -lt $3 ] || [ ${rightlen} -lt $3 ]; then   #filter out the kmer with at least one flank of <$flnk_len
grep ${mykmer}_ flank_coor.txt >> kmer_flanktooshort_flkcoor.txt  #store the flank coordinates of kmer with too short flank
echo ${mykmer} >> kmer_flanktooshort_4rm.txt   #storing the kmer with flank being too short for remove
fi

done  < header.txt #closing for looping through sig kmers header lines]

#remove the kmer with flanks being too short from the muktifasta file for blasting
awk 'BEGIN{while((getline<"kmer_flanktooshort_4rm.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' $1 > kmer_forblast.fasta 


###blasting: querys are kmers in kmer_forblast.fasta; subjects are genomes in gen_input.fasta, 
blastn -query kmer_forblast.fasta -subject $2 -outfmt 6 -out myout.txt
