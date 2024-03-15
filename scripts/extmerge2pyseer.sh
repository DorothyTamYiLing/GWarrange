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
                *)
                   echo "$1 is not a recognized flag!"
                   return 1;
                   ;;
          esac
  done  

echo "checking extmerge2pyseer_script parameters"
echo "extend_para: $extend_para";
echo "merge_para: $merge_para";
echo "klen: $klen";
echo "minor: $minor";
echo "major: $major";
echo "pheno: $pheno";
echo "pyseer_arg: $pyseer_arg";
echo "fsmlite_arg: $fsmlite_arg";
echo "unitigcaller_arg: $unitigcaller_arg";
echo "string_type: $string_type";

#make genome set using extension and merging parameters
echo "making genome set with extension and merge parameters"
Rscript scripts/merge_IS.R --input blastrep_out.txt --extend ${extend_para} --merge ${merge_para}
python3 scripts/iSreplace_2col.py --input fixed_genomes.fasta --coor ext${extend_para}_merge${merge_para}_mergedIS.txt --out ext${extend_para}_merge${merge_para}_ISreplaced_genomes

#enter into the directory of genome set
cd ext${extend_para}_merge${merge_para}_ISreplaced_genomes

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
${pyseer_arg} > kmer_pyseer"

#get p threshold
python3 ../scripts/count_patterns.py kmer_patterns.txt > count_pattern.txt

#extract threshold from pattern file
mythres=$(grep "Threshold" count_pattern.txt | cut -f2 -d$':' )

#get sigk
awk -v var="${mythres}" '{ if ($4 <= var) { print } }' kmer_pyseer > sigk_pyseer.txt #kmer

#convert to fasta format
awk '{print $1}' sigk_pyseer.txt > sigk_seq.txt 

#create multifasta file for significant kmer sequences
number=$(cat sigk_seq.txt | wc -l)

#remove any existing header file
if ! [ -f header.txt ]; then
  echo "header file does not exist."
else
  rm header.txt
fi

START=1
let "END=$number" 
 
for (( c=$START; c<=$END; c++ ))
do
    echo ">kmer""$c " >> header.txt
done

paste -d \\n header.txt sigk_seq.txt > final_sig.fasta

#classifying kmers into those containing "N" and those do not contain "N". Produce "sigk_noN.fasta" and "sigk_withN.fasta".
python3 ../scripts/class_k.py --input final_sig.fasta --outdir .

#rename
mv sigk_noN.fasta classk_noN.fasta 
mv sigk_withN.fasta classk_withN.fasta



#process unitigs
if [ ${string_type} = "kmer_and_unitigs" ]
then

echo "generating unitig-caller input file"
#generating unitig-caller input file
ls -d -1 $PWD/*_ISreplaced.fasta > input.txt
#running unitig-caller
eval "unitig-caller --call --pyseer --refs input.txt --out unitigcall_out ${unitigcaller_arg}"
#fixing sample names in output file for pyseer
sed "s/_ext${extend_para}_merge${merge_para}_ISreplaced//g" unitigcall_out.pyseer > unitigcall.pyseer
gzip unitigcall.pyseer

#run pyseer on unitig-caller output
echo "run pyseer on unitig-caller output"
eval "pyseer --phenotypes ${pheno} \
--kmers unitigcall.pyseer.gz \
--output-patterns unitig_patterns.txt \
${pyseer_arg} > unitig_pyseer"

#get p threshold
python3 ../scripts/count_patterns.py unitig_patterns.txt > count_uni_pattern.txt

#extract threshold from pattern file
mythres_uni=$(grep "Threshold" count_uni_pattern.txt | cut -f2 -d$':' )

#get sig unitig
awk -v var="${mythres_uni}" '{ if ($4 <= var) { print } }' unitig_pyseer > siguni_pyseer.txt #unitig

#convert to fasta format, unitigs
awk '{print $1}' siguni_pyseer.txt > siguni_seq.txt 

#create multifasta file for significant unitig sequences
number=$(cat siguni_seq.txt | wc -l)

#remove any existing header file
if ! [ -f header.txt ]; then
  echo "header file does not exist."
else
  rm header.txt
fi

START=1
let "END=$number" 
 
for (( c=$START; c<=$END; c++ ))
do
    echo ">unitig""$c " >> header.txt
done

paste -d \\n header.txt siguni_seq.txt > siguni_seq.fasta

#add kmer names to pyseer output
sed 's/>//g' header.txt > header_1.txt
paste header_1.txt sigk_pyseer.txt > sigk_pyseer_1.txt & mv sigk_pyseer_1.txt sigk_pyseer.txt
rm header_1.txt

#add header line
pyhead=$(head -1 kmer_pyseer)
pyhead_1="kmer_ID   "${pyhead}   
echo ${pyhead_1} >  pyhead.txt
cat pyhead.txt sigk_pyseer.txt > sigk_pyseer_1.txt mv sigk_pyseer_1.txt sigk_pyseer.txt

#classifying unitigs into those containing "N" and those do not contain "N". Produce "sigk_noN.fasta" and "sigk_withN.fasta".
python3 ../scripts/class_k.py --input siguni_seq.fasta --outdir .

# concatenatingn sig. kmers with N and sig. unitigs, final_sig.fasta will overwrite the previous one from kmers
cat classk_withN.fasta siguni_seq.fasta > final_sig.fasta

fi


cd ..
