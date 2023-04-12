#set default parameters in getopts
flk_len=30
dedupk=2
intkrd=1000

while getopts k:g:p:l:f:d:o:s:d:r: flag
do
    case "${flag}" in
        k) sigk=${OPTARG};;
        g) gen=${OPTARG};;
        p) pheno=${OPTARG};;
        l) k_len=${OPTARG};;
        f) flk_len=${OPTARG};;
        d) flk_dist=${OPTARG};;
        o) outdir=${OPTARG};;
        s) gen_size=${OPTARG};;
        x) dedupk==${OPTARG};;
        y) intkrd==${OPTARG};;
    esac
done
echo "sigk: $sigk";
echo "gen: $gen";
echo "pheno: $pheno";
echo "k_len: $k_len";
echo "flk_len: $flk_len";
echo "flk_dist: $flk_dist";
echo "outdir: $outdir";
echo "genome size for plot: $gen_size";
echo "genome pos sig digit for dedup split k in plot: $dedupk";
echo "round off intact k genome position to the nearest multiple of an integar: $intkrd";

#create output directory
mkdir $outdir

#python3 class_k.py --input ./example_data/clus1clus2_sigk_withN.fasta --outdir ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir
python3 ./scripts/class_k.py --input $sigk --outdir $outdir
echo "classifying sig k"

gunzip $gen
echo "unzipping genomes"

#processing sig kmers with N, only when the input fasta file is present
[ -s ${outdir}/sigk_withN.fasta ] && withN=1 || withN=0

if [[ ${withN} -eq 1 ]]
then
echo "there is sigk withN"
bash ./scripts/filtering_kmer_and_blast.sh ${outdir}/sigk_withN.fasta ${gen/.gz} $outdir $flk_len
echo "now run make_flank_summary.R"
Rscript ./scripts/make_flank_summary.R --k.len $k_len --pheno $pheno --outdir $outdir --flkdist $flk_dist  --dedupk ${dedupk}
else
echo "no withN sig k"
fi

#processing sig kmers without N, only when the input fasta file is present
[ -s ${outdir}/sigk_noN.fasta ] && noN=1 || noN=0

if [[ ${noN} -eq 1 ]]
then
echo  "there is sigk_noN.fasta"
blastn -query ${outdir}/sigk_noN.fasta -subject ${gen/.gz} -outfmt 6 -out ${outdir}/mynoN_out.txt
echo "now run process_sigkNoN.R"
Rscript ./scripts/process_sigkNoN.R --pheno $pheno --outdir $outdir
else
echo "no NoN sig k"
fi

#plot split kmers
echo "plotting split kmers"
mysplitkdir=${outdir}/splitk_plots
mkdir ${mysplitkdir}
sed 1d ${outdir}/myshort_splitk_out_uniq.txt > ${outdir}/myshort_splitk_out_uniq_nohead.txt
for x in $(cut -f1 ${outdir}/myshort_splitk_out_uniq_nohead.txt)
do
 echo $x
 Rscript ./scripts/plot_flk_kmer_prop.R --kmer $x --phen $pheno --coor ${outdir}/myflk_behave_pheno.txt --genome.size ${gen_size} --outdir ${mysplitkdir}/plot${x} --flk.dist ${flk_dist}
done
rm ${outdir}/myshort_splitk_out_uniq_nohead.txt

#plot intact k
[ -s ${outdir}/myintactkwithN_out.txt ] && intkN=1 || intkN=0
if [[ ${intkN} -eq 1 ]]
then

echo "plotting intact k with N"

#subset the intact kmers that are reverse in majority of the 0 genomes, and forward and majority of the 1 genomes
awk -F "\t" 'NR==1; NR > 1{ if ($6 < 0.5 && $7 > 0.5 && $8 > 0.5 && $9 < 0.5) { print } }' ${outdir}/myintactkwithN_out.txt > ${outdir}/myintactkwithN_rev0fwd1_set.txt

#subset the intact kmers that are reverse in majority of the 1 genomes, and forward and majority of the 0 genomes
awk -F "\t" 'NR==1; NR > 1{ if ($6 > 0.5 && $7 < 0.5 && $8 < 0.5 && $9 > 0.5) { print } }' ${outdir}/myintactkwithN_out.txt > ${outdir}/myintactkwithN_rev1fwd0_set.txt

Rscript ./scripts/plot_intactk.R --input ${outdir}/myintactkwithN_rev0fwd1_set.txt \
--outdir ${outdir} \
--outname myintactkwithN_rev0fwd1 \
--gen_size ${gen_size} \
--intkrd ${intkrd}

Rscript ./scripts/plot_intactk.R --input ${outdir}/myintactkwithN_rev1fwd0_set.txt \
--outdir ${outdir} \
--outname myintactkwithN_rev1fwd0 \
--gen_size ${gen_size} \
--intkrd ${intkrd}

else
echo "no myintactkwithN_out.txt for plot"
fi

[ -s ${outdir}/myNoNintactk_out.txt ] && intknoN=1 || intknoN=0
if [[ ${intknoN} -eq 1 ]]
then

echo "plotting intact k withoutN"

#subset the intact kmers that are reverse in majority of the 0 genomes, and forward and majority of the 1 genomes
awk -F "\t" 'NR==1; NR > 1{ if ($6 < 0.5 && $7 > 0.5 && $8 > 0.5 && $9 < 0.5) { print } }' ${outdir}/myNoNintactk_out.txt > ${outdir}/myNoNintactk_rev0fwd1_set.txt

#subset the intact kmers that are reverse in majority of the 1 genomes, and forward and majority of the 0 genomes
awk -F "\t" 'NR==1; NR > 1{ if ($6 > 0.5 && $7 < 0.5 && $8 < 0.5 && $9 > 0.5) { print } }' ${outdir}/myNoNintactk_out.txt > ${outdir}/myNoNintactk_rev1fwd0_set.txt

Rscript ./scripts/plot_intactk.R --input ${outdir}/myNoNintactk_rev0fwd1_set.txt \
--outdir ${outdir} \
--outname myNoNintactk_rev0fwd1 \
--gen_size ${gen_size} \
--intkrd ${intkrd}

Rscript ./scripts/plot_intactk.R --input ${outdir}/myNoNintactk_rev1fwd0_set.txt \
--outdir ${outdir} \
--outname myNoNintactk_rev1fwd0 \
--gen_size ${gen_size} \
--intkrd ${intkrd}

else
echo "no myintactknoN_out.txt for plot"
fi



gzip ${gen/.gz}

