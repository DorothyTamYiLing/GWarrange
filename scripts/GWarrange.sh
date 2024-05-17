#set default parameters in getopts
flk_len=30
dedupk=2
intkrd=1000
ext_mrg_min="100_3"
ext_mrg_max="7000_3"
pyseer_arg="--min-af 0.05 --max-af 0.95"
fsmlite_arg="-v -t tmp -m 200 -M 200"
unitigcaller_arg=""
string_type="kmer"
thread=8
exp_fac=86
yaxis=360
arr_dist=70
split_h=7
split_w=10
merge=40000
intact_h=100
intact_w=180
intact_res=150

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
                -thread)
                    shift
                    thread=$1
                    shift
                    ;;
                -exp_fac)
                    shift
                    exp_fac=$1
                    shift
                    ;;
                -yaxis)
                    shift
                    yaxis=$1
                    shift
                    ;;
                -arr_dist)
                    shift
                    arr_dist=$1
                    shift
                    ;;
                -split_h)
                    shift
                    split_h=$1
                    shift
                    ;;
                -split_w)
                    shift
                    split_w=$1
                    shift
                    ;;
                -merge)
                    shift
                    merge=$1
                    shift
                    ;;
                -intact_h)
                    shift
                    intact_h=$1
                    shift
                    ;;
                -intact_w)
                    shift
                    intact_w=$1
                    shift
                    ;;
                -intact_res)
                    shift
                    intact_res=$1
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
echo "number of thread for BLAST, pyseer and unitig-callers: $thread";
echo "minimum extension and merge parameters: $ext_mrg_min";
echo "maximum extension and merge parameters: $ext_mrg_max";
echo "how much the arrow expand horizontally for visibility in relative to the genome size: $exp_fac";
echo "y-axis height of the plot: $yaxis";
echo "vertical distance between arrows: $arr_dist";
echo "height of the split-k plot: $split_h";
echo "width of the split-k plot: $split_w";
echo "merge arrows into one when they are less than XXXbp apart: $merge";
echo "height of the intact-k plot: $intact_h";
echo "width of the intact-k plot: $intact_w";
echo "set intact k plot resolution: $intact_res";


min_ext=$(echo ${ext_mrg_min} | cut -d "_" -f1)
min_mrg=$(echo ${ext_mrg_min} | cut -d "_" -f2)
max_ext=$(echo ${ext_mrg_max} | cut -d "_" -f1)
max_mrg=$(echo ${ext_mrg_max} | cut -d "_" -f2)
echo "min_ext: $min_ext"
echo "min_mrg: $min_mrg"
echo "max_ext: $max_ext"
echo "max_mrg: $max_mrg"

#adding header to phenotype file for pyseer input format
echo "samples binary" | cat - $pheno > ${pheno/.txt}_4pyseer.txt

#define genome file name
genname=$(basename ${gen/.fna.gz})
echo ${genname}

#unzip the genome file if neccesasry

gunzip $gen

echo "blast starting gene with genomes"

#blast start gene with genomes

blastn -query ${startgene} \
-subject ${gen/.gz} \
-outfmt 6 -out blast_startgene_out.txt

echo "reorientating genomes with the starting gene"

#reorientating genomes with the starting gene
python3 scripts/fix_genome.py --input ${gen/.gz} --mycoor blast_startgene_out.txt

echo "get coordinates for repeat sequences in reorientated genomes"

#get coordinates for repeat sequences in reorientated genomes
blastn -query $replist \
-subject fixed_genomes.fasta \
-outfmt 6 -out blastrep_out.txt

############ script for running extension+merge, fsm-lite, unitig-caller, pyseer, pyseer post-processing ########################

#build database for efficient blasting
echo "build database for efficient blasting"
makeblastdb -in ${gen/.gz} -dbtype nucl -out genome_db

#if both max and min parameter are the same, then just perform one set
if [ "${ext_mrg_min}" != "${ext_mrg_max}" ]

then

echo "run with both minimum and maximum extension and merge parameters"

echo "Using minimum extension and merge parameters"

bash scripts/extmerge2pyseer.sh \
-extend_para ${min_ext} -merge_para ${min_mrg} \
-pheno ${pheno/.txt}_4pyseer.txt \
-pyseer_arg "${pyseer_arg}" \
-fsmlite_arg "${fsmlite_arg}" \
-unitigcaller_arg "${unitigcaller_arg}" \
-string_type ${string_type} \
-thread ${thread}

#calculate value for d flag
min_cut=$(grep "Maximum size" ext${min_ext}_merge${min_mrg}_mergedISstat.txt | cut -d$'\t' -f 2)
echo "Maximum size of merged repeats:"${min_cut}
min_d=$(awk "BEGIN { print int(${min_cut}*1.2)}")
echo "value for -d flag (1.2 times of above value) :"${min_d}

if [ -s ext${min_ext}_merge${min_mrg}_ISreplaced_genomes_${string_type}/final_sig.fasta ]; then
    echo "process final significant kmers/unitigs for detecting rearrangements"
#running main.sh For ext100_merge3_ISreplaced_genomes set unitigs
bash scripts/main.sh -sigk ext${min_ext}_merge${min_mrg}_ISreplaced_genomes_${string_type}/final_sig.fasta \
-gen genome_db \
-pheno ${pheno} -flk_dist ${min_d} -flk_len ${flk_len} \
-outdir ${genname}_ext${min_ext}_merge${min_mrg}_${string_type}_outdir \
-gen_size ${gen_size} -dedupk ${dedupk} -intkrd ${intkrd} -thread ${thread} \
-exp_fac ${exp_fac} -yaxis ${yaxis} -arr_dist ${arr_dist} -split_h ${split_h} -split_w ${split_w} -merge ${merge} \
-intact_h ${intact_h} -intact_w ${intact_w} -intact_res ${intact_res}
else
    echo "no significant kmer/unitig"
fi

echo "Using maximum extension and merge parameters"

bash scripts/extmerge2pyseer.sh \
-extend_para ${max_ext} -merge_para ${max_mrg} \
-pheno ${pheno/.txt}_4pyseer.txt \
-pyseer_arg "${pyseer_arg}" \
-fsmlite_arg "${fsmlite_arg}" \
-unitigcaller_arg "${unitigcaller_arg}" \
-string_type ${string_type} \
-thread ${thread}


#calculate value for d flag
max_cut=$(grep "Maximum size" ext${max_ext}_merge${max_mrg}_mergedISstat.txt | cut -d$'\t' -f 2)
echo "Maximum size of merged repeats:"${max_cut}
min_d=$(awk "BEGIN { print int(${max_cut}*1.2)}")
echo "value for -d flag (1.2 times of above value) :"${min_d}

if [ -s ext${max_ext}_merge${max_mrg}_ISreplaced_genomes_${string_type}/final_sig.fasta ]; then
    echo "process final significant kmers/unitigs for detecting rearrangements"
#running main.sh For ext7000_merge200_ISreplaced_genomes set unitigs
bash scripts/main.sh -sigk ext${max_ext}_merge${max_mrg}_ISreplaced_genomes_${string_type}/final_sig.fasta \
-gen genome_db \
-pheno ${pheno} -flk_dist ${min_d} -flk_len ${flk_len} \
-outdir ${genname}_ext${max_ext}_merge${max_mrg}_${string_type}_outdir \
-gen_size ${gen_size} -dedupk ${dedupk} -intkrd ${intkrd} -thread ${thread} \
-exp_fac ${exp_fac} -yaxis ${yaxis} -arr_dist ${arr_dist} -split_h ${split_h} -split_w ${split_w} -merge ${merge} \
-intact_h ${intact_h} -intact_w ${intact_w} -intact_res ${intact_res}
else
     echo "no significant kmer/unitig"
fi

else  # else for if [ "${ext_mrg_min}" != "${ext_mrg_max}" ]

echo "run with one set of extension and merge parameters"

bash scripts/extmerge2pyseer.sh \
-extend_para ${min_ext} -merge_para ${min_mrg} \
-pheno ${pheno/.txt}_4pyseer.txt \
-pyseer_arg "${pyseer_arg}" \
-fsmlite_arg "${fsmlite_arg}" \
-unitigcaller_arg "${unitigcaller_arg}" \
-string_type ${string_type} \
-thread ${thread}

#calculate value for d flag
min_cut=$(grep "Maximum size" ext${min_ext}_merge${min_mrg}_mergedISstat.txt | cut -d$'\t' -f 2)
echo "Maximum size of merged repeats:"${min_cut}
min_d=$(awk "BEGIN { print int(${min_cut}*1.2)}")
echo "value for -d flag (1.2 times of above value) :"${min_d}

if [ -s ext${min_ext}_merge${min_mrg}_ISreplaced_genomes_${string_type}/final_sig.fasta ]; then
    echo "process final significant kmers/unitigs for detecting rearrangements"
#running main.sh For ext100_merge3_ISreplaced_genomes set unitigs
bash scripts/main.sh -sigk ext${min_ext}_merge${min_mrg}_ISreplaced_genomes_${string_type}/final_sig.fasta \
-gen genome_db \
-pheno ${pheno} -flk_dist ${min_d} -flk_len ${flk_len} \
-outdir ${genname}_ext${min_ext}_merge${min_mrg}_${string_type}_outdir \
-gen_size ${gen_size} -dedupk ${dedupk} -intkrd ${intkrd} -thread ${thread} \
-exp_fac ${exp_fac} -yaxis ${yaxis} -arr_dist ${arr_dist} -split_h ${split_h} -split_w ${split_w} -merge ${merge} \
-intact_h ${intact_h} -intact_w ${intact_w} -intact_res ${intact_res}
else
     echo "no significant kmer/unitig"
fi

fi   #close bracket for if [ "${ext_mrg_min}" != "${ext_mrg_max}" 

echo "DONE"