#set default parameters in getopts
flk_len=30
dedupk=2
intkrd=1000
ext_mrg_min="100_3"
ext_mrg_max="7000_200"
pyseer_arg="--min-af 0.05 --max-af 0.95 --print-samples --no-distances --cpu 8"
fsmlite_arg="-v -t tmp -m 200 -M 200"
unitigcaller_arg="--threads 8"
string_type="kmer_and_unitigs"
thread_blast=8

while test $# -gt 0; do
           case "$1" in
                -gen)
                    shift
                    gen=$1
                    shift
                    ;;
                -pheno)
                    shift
                    pheno=$1
                    shift
                    ;;
                -flk_len)
                    shift
                    flk_len=$1
                    shift
                    ;; 
                -gen_size)
                    shift
                    gen_size=$1
                    shift
                    ;;
                -dedupk)
                    shift
                    dedupk=$1
                    shift
                    ;;
                -intkrd)
                    shift
                    intkrd=$1
                    shift
                    ;;
                -startgene)
                    shift
                    startgene=$1
                    shift
                    ;;
                -replist)
                    shift
                    replist=$1
                    shift
                    ;;
                -ext_mrg_min)
                    shift
                    ext_mrg_min=$1
                    shift
                    ;;
		-ext_mrg_max)
                    shift
                    ext_mrg_max=$1
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
                    uniticaller_arg=$1
                    shift
                    ;; 
                -string_type)
                    shift
                    string_type=$1
                    shift
                    ;;
                -thread_blast)
                    shift
                    thread_blast=$1
                    shift
                    ;;
                *)
                   echo "$1 is not a recognized flag!"
                   return 1;
                   ;;
          esac
  done  

echo "gen: $gen";
echo "pheno: $pheno";
echo "flk_len: $flk_len";
echo "genome size for plot: $gen_size";
echo "genome pos sig digit for dedup split k in plot: $dedupk";
echo "round off intact k genome position to the nearest multiple of an integar: $intkrd";
echo "starting gene: $startgene";
echo "repeat sequence list: $replist";
echo "pyseer_arg: $pyseer_arg";
echo "fsmlite_arg: $fsmlite_arg";
echo "unitigcaller_arg: $unitigcaller_arg";
echo "string_type: $string_type";
echo "number of thread for BLAST: $thread_blast";
echo "minimum extension and merge parameters: $ext_mrg_min";
echo "maximum extension and merge parameters: $ext_mrg_max";

#ext_mrg="100_3"
#pheno="prn_status_pheno.txt"
#gen="PRN_468.fna.gz"
#startgene="gidA.fasta"
#replist="IS_NZ_CP025371.1.fasta"
#string_type="kmer"


#adding header to phenotype file for pyseer input format
echo "samples binary" | cat - example_data/$pheno > example_data/${pheno/.txt}_4pyseer.txt


#unzip the genome file if neccesasry
gunzip example_data/$gen

echo "blast gidA with genomes"

#blast gidA with genomes
blastn -query example_data/${startgene} \
-subject example_data/${gen/.gz} \
-outfmt 6 -out ${gen/.fna.gz}_${startgene/.fasta}_out.txt

echo "reorientating genomes with the starting gene"

#reorientating genomes with the starting gene
python3 scripts/fix_genome.py --input example_data/${gen/.gz} --mycoor ${gen/.fna.gz}_${startgene/.fasta}_out.txt

echo "get coordinates for repeat sequences in reorientated genomes"

#get coordinates for repeat sequences in reorientated genomes
blastn -query example_data/$replist \
-subject fixed_genomes.fasta \
-outfmt 6 -out blastrep_out.txt

############ script for running extension+merge, fsm-lite, unitig-caller, pyseer, pyseer post-processing ########################

echo "${pyseer_arg}"
echo "${fsmlite_arg}"
echo "${unitigcaller_arg}"

echo "run extmerge2pyseer.sh using minimum extend merge parameters"

min_ext=$(echo ${ext_mrg_min} | cut -d "_" -f1)
min_mrg=$(echo ${ext_mrg_min} | cut -d "_" -f2)
max_ext=$(echo ${ext_mrg_max} | cut -d "_" -f1)
max_mrg=$(echo ${ext_mrg_max} | cut -d "_" -f2)

bash scripts/extmerge2pyseer.sh \
-extend_para ${min_ext} -merge_para ${min_mrg} \
-pheno ../example_data/${pheno/.txt}_4pyseer.txt \
-pyseer_arg "${pyseer_arg}" \
-fsmlite_arg "${fsmlite_arg}" \
-unitigcaller_arg "${unitigcaller_arg}" \
-string_type ${string_type}

echo "run extmerge2pyseer.sh using maximum extend merge parameters"

bash scripts/extmerge2pyseer.sh \
-extend_para ${max_ext} -merge_para ${max_mrg} \
-pheno ../example_data/${pheno/.txt}_4pyseer.txt \
-pyseer_arg "${pyseer_arg}" \
-fsmlite_arg "${fsmlite_arg}" \
-unitigcaller_arg "${unitigcaller_arg}" \
-string_type ${string_type}


echo "build database for efficient blasting"
#build database for efficient blasting
makeblastdb -in example_data/${gen/.gz} -dbtype nucl -out genome_db

#calculate value for d flag
min_cut=$(grep "Maximum size" ext${min_ext}_merge${min_mrg}_mergedISstat.txt | cut -d$'\t' -f 2)
min_d=$(awk "BEGIN { print int(${min_cut}*1.2)}")

max_cut=$(grep "Maximum size" ext${max_ext}_merge${max_mrg}_mergedISstat.txt | cut -d$'\t' -f 2)
max_d=$(awk "BEGIN { print int(${max_cut}*1.2)}")


#running main.sh For ext100_merge3_ISreplaced_genomes set unitigs
bash scripts/main.sh -k ext${min_ext}_merge${min_mrg}_ISreplaced_genomes/final_sig.fasta \
-g genome_db \
-p example_data/${pheno} -d ${min_d} -f ${flk_len} \
-o ${gen/.fna.gz}_ext${min_ext}_merge${min_mrg}_outdir -s ${gen_size} -x ${dedupk} -y ${intkrd} -t ${thread_blast}

#running main.sh For ext7000_merge200_ISreplaced_genomes set unitigs
bash scripts/main.sh -k ext${max_ext}_merge${max_mrg}_ISreplaced_genomes/final_sig.fasta \
-g genome_db \
-p example_data/${pheno} -d ${max_d} -f ${flk_len} \
-o ${gen/.fna.gz}_ext${max_ext}_merge${max_mrg}_outdir -s ${gen_size} -x ${dedupk} -y ${intkrd} -t ${thread_blast}


