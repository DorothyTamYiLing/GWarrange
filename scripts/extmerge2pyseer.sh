# extmerge2pyseer.sh #
while test $# -gt 0; do
           case "$1" in
                -extend_para)
                    shift
                    extend_para=$1
                    shift
                    ;;
                -merge_para)
                    shift
                    merge_para=$1
                    shift
                    ;;
                -klen)
                    shift
                    klen=$1
                    shift
                    ;;
                -minor)
                    shift
                    minor=$1
                    shift
                    ;;
                -major)
                    shift
                    major=$1
                    shift
                    ;;
                -pheno)
                    shift
                    pheno=$1
                    shift
                    ;;
                -pyseer_arg)
                    shift
                    pyseer_arg=$1
                    shift
                    ;;
		-fsmlite_arg)
                    shift
                    fsmlite_arg=$1
                    shift
                    ;;
		-unitigcaller_arg)
                    shift
                    unitigcaller_arg=$1
                    shift
                    ;;	 
               -string_type)
                    shift
                    string_type=$1
                    shift
		    ;;
            -thread)
                    shift
                    thread=$1
                    shift
            ;;
                *)
                   echo "$1 is not a recognized flag!"
                   return 1;
                   ;;
          esac
  done  

#echo "checking extmerge2pyseer_script parameters"
echo "extend_para: $extend_para";
echo "merge_para: $merge_para";
#echo "pyseer_arg: $pyseer_arg";
#echo "fsmlite_arg: $fsmlite_arg";
#echo "unitigcaller_arg: $unitigcaller_arg";
#echo "string_type: $string_type";
#echo "thread: $thread";

#make genome set using extension and merging parameters
echo "making genome set with extension and merge parameters"
Rscript scripts/merge_IS.R --input blastrep_out.txt --extend ${extend_para} --merge ${merge_para}

#create output directory, replace old one if exists
if [[ -d ext${extend_para}_merge${merge_para}_ISreplaced_genomes_${string_type} ]]; then
        echo "directory exists, replacing with the new one"
        rm -r ext${extend_para}_merge${merge_para}_ISreplaced_genomes_${string_type}
        python3 scripts/iSreplace_2col.py --input fixed_genomes.fasta --coor ext${extend_para}_merge${merge_para}_mergedIS.txt --out ext${extend_para}_merge${merge_para}_ISreplaced_genomes
else
        python3 scripts/iSreplace_2col.py --input fixed_genomes.fasta --coor ext${extend_para}_merge${merge_para}_mergedIS.txt --out ext${extend_para}_merge${merge_para}_ISreplaced_genomes

fi


#rename to add string type
mv ext${extend_para}_merge${merge_para}_ISreplaced_genomes ext${extend_para}_merge${merge_para}_ISreplaced_genomes_${string_type}

#enter into the directory of genome set
cd ext${extend_para}_merge${merge_para}_ISreplaced_genomes_${string_type}

#process kmers
#generating fsm-ite input file
echo "generating fsm-ite input file"
for f in *_ext${extend_para}_merge${merge_para}_ISreplaced.fasta; do id=$(basename "$f" _ext${extend_para}_merge${merge_para}_ISreplaced.fasta); echo $id $f; done > input.list
eval "fsm-lite -l input.list ${fsmlite_arg} | gzip - > fsm_output.txt.gz"

#run pyseer on fsm-lite output
echo "run pyseer on fsm-lite output"
eval "pyseer --phenotypes ${pheno} \
--kmers fsm_output.txt.gz \
--output-patterns kmer_patterns.txt \
${pyseer_arg} --cpu ${thread} > kmer_pyseer"

#get p threshold
python3 ../scripts/count_patterns.py kmer_patterns.txt > count_pattern.txt

#extract threshold from pattern file
mythres=$(grep "Threshold" count_pattern.txt | cut -f2 -d$':' )

#get sigk
awk -v var="${mythres}" '{ if ($4 <= var) { print } }' kmer_pyseer > sigk_pyseer.txt #kmer


#check if the file is empty
if [ -s sigk_pyseer.txt ]; then

sigkcount=$(wc -l sigk_pyseer.txt)

myk="yes"

echo "there are "${sigkcount}" significant kmers found"

#convert to fasta format
awk '{print $1}' sigk_pyseer.txt > sigk_seq.txt 

#create multifasta file for significant kmer sequences
number=$(cat sigk_seq.txt | wc -l)

#remove any existing header file
if ! [ -f header.txt ]; then
  echo "there is no header file"
else
    echo "remove existing header file"
  rm header.txt
fi

seq 1 ${number} | awk '$1=">kmer"$1' > header.txt

paste -d \\n header.txt sigk_seq.txt > sigk_seq.fasta

#create the final fasta, can overwrite it later
cp sigk_seq.fasta final_sig.fasta

#for combine with sig uni classifying kmers into those containing "N" and those do not contain "N". Produce "sigk_noN.fasta" and "sigk_withN.fasta"
#there is always output, just may produce empty file
python3 ../scripts/class_k.py --input sigk_seq.fasta --outdir .

sigkN=$(wc -l sigk_withN.fasta)
sigknoN=$(wc -l sigk_noN.fasta)

echo "there are "${sigkN}" significant kmers with N"
echo "there are "${sigknoN}" significant kmers without N"

#rename
mv sigk_noN.fasta classk_noN.fasta 
mv sigk_withN.fasta classk_withN.fasta

#add kmer names to pyseer output
sed 's/>//g' header.txt > header_1.txt
paste header_1.txt sigk_pyseer.txt > sigk_pyseer_1.txt
rm header_1.txt

#add header line
pyhead=$(head -1 kmer_pyseer)
pyhead_1="kmer_ID   "${pyhead}   
echo ${pyhead_1} >  pyhead.txt
cat pyhead.txt sigk_pyseer_1.txt > sigk_pyseer.txt
rm pyhead.txt
rm sigk_pyseer_1.txt


else
    echo "no significant kmer found"
    myk="no"
fi  #close bracket for if [ -s sigk_pyseer.txt ]; then



#################################

#process unitigs
if [ ${string_type} = "kmers_and_unitigs" ]
then

echo "generating unitig-caller input file"
#generating unitig-caller input file
ls -d -1 $PWD/*_ISreplaced.fasta > input.txt
#running unitig-caller
eval "unitig-caller --call --pyseer --refs input.txt --out unitigcall_out --threads ${thread} ${unitigcaller_arg}"
#fixing sample names in output file for pyseer
sed "s/_ext${extend_para}_merge${merge_para}_ISreplaced//g" unitigcall_out.pyseer > unitigcall.pyseer
gzip unitigcall.pyseer

#run pyseer on unitig-caller output
echo "run pyseer on unitig-caller output"
eval "pyseer --phenotypes ${pheno} \
--kmers unitigcall.pyseer.gz \
--output-patterns unitig_patterns.txt \
${pyseer_arg} --cpu ${thread} > unitig_pyseer"

#get p threshold
python3 ../scripts/count_patterns.py unitig_patterns.txt > count_uni_pattern.txt

#extract threshold from pattern file
mythres_uni=$(grep "Threshold" count_uni_pattern.txt | cut -f2 -d$':' )

#get sig unitig
awk -v var="${mythres_uni}" '{ if ($4 <= var) { print } }' unitig_pyseer > siguni_pyseer.txt #unitig

#check if the file is filled
if [ -s siguni_pyseer.txt ]; then

myuni="yes"

sigunicount=$(wc -l siguni_pyseer.txt)

echo "there are "${sigunicount}" significant unitigs found"

#convert to fasta format, unitigs
awk '{print $1}' siguni_pyseer.txt > siguni_seq.txt 

#create multifasta file for significant unitig sequences
number=$(cat siguni_seq.txt | wc -l)

#remove any existing header file
if ! [ -f header.txt ]; then
  echo "there is no header file"
else
    echo "remove existing header file"
  rm header.txt
fi

seq 1 ${number} | awk '$1=">unitig"$1' > header.txt

paste -d \\n header.txt siguni_seq.txt > siguni_seq.fasta

#classifying unitigs into those containing "N" and those do not contain "N". Produce "sigk_noN.fasta" and "sigk_withN.fasta".
#python3 ../scripts/class_k.py --input siguni_seq.fasta --outdir .

#add unitig names to pyseer output
sed 's/>//g' header.txt > header_1.txt
paste header_1.txt siguni_pyseer.txt > siguni_pyseer_1.txt
rm header_1.txt

#add header line
pyhead=$(head -1 unitig_pyseer)
pyhead_1="kmer_ID   "${pyhead}   
echo ${pyhead_1} >  pyhead.txt
cat pyhead.txt siguni_pyseer_1.txt > siguni_pyseer.txt
rm pyhead.txt
rm siguni_pyseer_1.txt

else

    echo "no significant unitig found"
    myuni="no"

fi #close bracket for if [ -s siguni_pyseer.txt ]; then


#echo $myk $myuni

#if there is sig kmer but no sig uni, do nothing, keeping using final_seq.fasta generated from before

#if there is sig kmer and sig uni, combine sigk with N and all sig uni, overwrite final_seq.fasta
if [ ${myk} = "yes" ] && [ ${myuni} = "yes" ]; then
    echo "there are both significant kmers and significant unitigs"
    echo "combining significant kmers with N and significant unitigs for rearrangement detection"
cat classk_withN.fasta siguni_seq.fasta > final_sig.fasta
fi

if [ ${myk} == "no" ] && [ ${myuni} == "no" ]; then
echo "no significant kmer nor significant unitig found"
touch final_sig.fasta  #create an empty file anyway
fi


fi #close bracket for if [ ${string_type} = "kmer_and_unitigs" ]

#reset everything
myk="none"
myuni="none"

cd ..

