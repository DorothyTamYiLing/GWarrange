#usage:bash filtering_kmer_and_blast.sh allsig_kmer_withN.fasta 111_yearGWAS_genlist.fasta.gz output 30

###filtering sig. kmers for blasting by keeping only the kmers that contain flanking sequences (both side) of at least $flnk_len bp in size

#get the flank start and end coordinates of the sig kmers, output file name: flank_coor.txt
python3 ./scripts/extract_flank_coor.py --input $1 --outdir $3

grep ">" $1 > $3/header.txt

if [[ -d $3/kmer_flanktooshort_flkcoor.txt ]]; then rm $3/kmer_flanktooshort_flkcoor.txt; fi

if [[ -d $3/kmer_flanktooshort_4rm.txt ]]; then rm $3/kmer_flanktooshort_4rm.txt; fi

#rm $3/kmer_flanktooshort_flkcoor.txt
#rm $3/kmer_flanktooshort_4rm.txt

#remove kmers with short flanks and extract kmers for blasting, by looping through the header file
while read mykmer
do
mykmer=$(echo ${mykmer} | sed 's/>//')  #remove ">" from the header line
#mykmer="kmer1" # specify which kmer being process
#echo ${mykmer}

leftflankend=$(grep ${mykmer}_ $3/flank_coor.txt | cut -d "_" -f2)
rightflankstart=$(grep ${mykmer}_ $3/flank_coor.txt | cut -d "_" -f3)
kmerlen=$(grep ${mykmer}_ $3/flank_coor.txt | cut -d "_" -f4)
rightlen="$((${kmerlen}-${rightflankstart}+1))"       

if [ ${leftflankend} -lt $4 ] || [ ${rightlen} -lt $4 ]; then   #filter out the kmer with at least one flank of <$flnk_len
grep ${mykmer}_ $3/flank_coor.txt >> $3/kmer_flanktooshort_flkcoor.txt  #store the flank coordinates of kmer with too short flank
echo ${mykmer} >> $3/kmer_flanktooshort_4rm.txt   #storing the kmer with flank being too short for remove
fi

done  < $3/header.txt #closing for looping through sig kmers header lines]

#remove the kmer with flanks being too short from the muktifasta file for blasting
name=$3/kmer_flanktooshort_4rm.txt
awk -v var="$name" 'BEGIN{while((getline<var)>0)l[">"$1]=1}/^>/{f=!l[$1]}f' $1 > $3/kmer_forblast.fasta 

#remove header.txt
rm $3/header.txt

###blasting: querys are kmers in kmer_forblast.fasta; subjects are genomes in gen_input.fasta, 
mygenome=$2
#gunzip  $mygenome
blastn -query $3/kmer_forblast.fasta -subject ${mygenome/.gz} -outfmt 6 -out $3/myout.txt
#gzip  ${mygenome/.gz}
