#set default parameters in getopts
flk_len=30
dedupk=2
intkrd=1000

while getopts k:g:p:f:d:o:s:d:x:y: flag
do
    case "${flag}" in
        k) sigk=${OPTARG};;
        g) gen=${OPTARG};;
        p) pheno=${OPTARG};;
        f) flk_len=${OPTARG};;
        d) flk_dist=${OPTARG};;
        o) outdir=${OPTARG};;
        s) gen_size=${OPTARG};;
        x) dedupk=${OPTARG};;
        y) intkrd=${OPTARG};;
    esac
done
echo "sigk: $sigk";
echo "gen: $gen";
echo "pheno: $pheno";
echo "flk_len: $flk_len";
echo "flk_dist: $flk_dist";
echo "outdir: $outdir";
echo "genome size for plot: $gen_size";
echo "genome pos sig digit for dedup split k in plot: $dedupk";
echo "round off intact k genome position to the nearest multiple of an integar: $intkrd";

#create output directory, replace old one if exists
if [[ -d ./$outdir ]]; then
        echo "directory exists, replacing with the new one"; rm -r ./$outdir; mkdir $outdir
else
        mkdir $outdir
fi

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
echo "now run filtering_kmer_and_blast.sh"
bash ./scripts/filtering_kmer_and_blast.sh ${outdir}/sigk_withN.fasta ${gen/.gz} $outdir $flk_len
echo "now run make_flank_summary.R"
Rscript ./scripts/make_flank_summary.R --pheno $pheno --outdir $outdir --flkdist $flk_dist  --dedupk ${dedupk}
else
echo "no withN sig k"
fi

#processing sig kmers without N, only when the input fasta file is present
[ -s ${outdir}/sigk_noN.fasta ] && noN=1 || noN=0

if [[ ${noN} -eq 1 ]]
then
echo  "there is sigk_noN.fasta"
python3 ./scripts/extract_knoN_length.py --input ${outdir}/sigk_noN.fasta --outdir ${outdir}
blastn -query ${outdir}/sigk_noN.fasta -subject ${gen/.gz} -outfmt 6 -out ${outdir}/mynoN_out.txt
echo "now run process_sigkNoN.R"
Rscript ./scripts/process_sigkNoN.R --pheno $pheno --outdir $outdir
else
echo "no NoN sig k"
fi

#plot split kmers
[ -s ${outdir}/myshort_splitk_out_uniq.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
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
else
echo "no split kmers for plot"
fi

#plot intact k
[ -s ${outdir}/myintactkwithN_out.txt ] && intkN=1 || intkN=0
if [[ ${intkN} -eq 1 ]]
then

echo "plotting intact k with N"
[ -s ${outdir}/myintactkwithN_rev0fwd1_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
echo "plotting myintactkwithN_rev0fwd1_set"
Rscript ./scripts/plot_intactk.R --input ${outdir}/myintactkwithN_rev0fwd1_set.txt \
--outdir ${outdir} \
--outname myintactkwithN_rev0fwd1 \
--gen_size ${gen_size} \
--intkrd ${intkrd}
fi

[ -s ${outdir}/myintactkwithN_rev1fwd0_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
echo "plotting myintactkwithN_rev1fwd0_set"
Rscript ./scripts/plot_intactk.R --input ${outdir}/myintactkwithN_rev1fwd0_set.txt \
--outdir ${outdir} \
--outname myintactkwithN_rev1fwd0 \
--gen_size ${gen_size} \
--intkrd ${intkrd}
fi

[ -s ${outdir}/myintactkwithN_other_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
echo "plotting myintactkwithN_other_set"
Rscript ./scripts/plot_intactk.R --input ${outdir}/myintactkwithN_other_set.txt \
--outdir ${outdir} \
--outname myintactkwithN_other \
--gen_size ${gen_size} \
--intkrd ${intkrd}
fi

else
echo "no myintactkwithN_out.txt for plot"
fi

[ -s ${outdir}/myNoNintactk_out.txt ] && intknoN=1 || intknoN=0
if [[ ${intknoN} -eq 1 ]]
then

echo "plotting intact k withoutN"

echo "plotting myNoNintactk_rev0fwd1"
[ -s ${outdir}/myNoNintactk_rev0fwd1_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
Rscript ./scripts/plot_intactk.R --input ${outdir}/myNoNintactk_rev0fwd1_set.txt \
--outdir ${outdir} \
--outname myNoNintactk_rev0fwd1 \
--gen_size ${gen_size} \
--intkrd ${intkrd}
fi

echo "plotting myNoNintactk_rev1fwd0"
[ -s ${outdir}/myNoNintactk_rev1fwd0_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
Rscript ./scripts/plot_intactk.R --input ${outdir}/myNoNintactk_rev1fwd0_set.txt \
--outdir ${outdir} \
--outname myNoNintactk_rev1fwd0 \
--gen_size ${gen_size} \
--intkrd ${intkrd}
fi

echo "plotting myNoNintactk_other"
[ -s ${outdir}/myNoNintactk_other_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
Rscript ./scripts/plot_intactk.R --input ${outdir}/myNoNintactk_other_set.txt \
--outdir ${outdir} \
--outname myNoNintactk_other \
--gen_size ${gen_size} \
--intkrd ${intkrd}
fi

else
echo "no myintactknoN_out.txt for plot"
fi

#plotting rev0fwd1 all intactk set
if [[ -f ${outdir}/myintactkwithN_rev0fwd1_set.txt && -f ${outdir}/myNoNintactk_rev0fwd1_set.txt ]]; then
sed '1d' ${outdir}/myNoNintactk_rev0fwd1_set.txt > ${outdir}/myNoNintactk_rev0fwd1_set_no1stline.txt
cat ${outdir}/myintactkwithN_rev0fwd1_set.txt ${outdir}/myNoNintactk_rev0fwd1_set_no1stline.txt > ${outdir}/myallintactk_rev0fwd1_set.txt
echo "plotting myallintactk_rev0fwd1_set"
Rscript ./scripts/plot_intactk.R --input ${outdir}/myallintactk_rev0fwd1_set.txt \
--outdir ${outdir} \
--outname myallintactk_rev0fwd1 \
--gen_size ${gen_size} \
--intkrd ${intkrd}
fi


#plotting rev1fwd0 all intactk set
if [[ -f ${outdir}/myintactkwithN_rev1fwd0_set.txt && -f ${outdir}/myNoNintactk_rev1fwd0_set.txt ]]; then
sed '1d' ${outdir}/myNoNintactk_rev1fwd0_set.txt > ${outdir}/myNoNintactk_rev1fwd0_set_no1stline.txt
cat ${outdir}/myintactkwithN_rev1fwd0_set.txt ${outdir}/myNoNintactk_rev1fwd0_set_no1stline.txt > ${outdir}/myallintactk_rev1fwd0_set.txt
echo "plotting myallintactk_rev1fwd0_set"
Rscript ./scripts/plot_intactk.R --input ${outdir}/myallintactk_rev1fwd0_set.txt \
--outdir ${outdir} \
--outname myallintactk_rev1fwd0 \
--gen_size ${gen_size} \
--intkrd ${intkrd}
fi

#placing output files into different output directories
mkdir ${outdir}/preprocssing

[ -s ${outdir}/flank_coor.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/flank_coor.txt ${outdir}/preprocssing
fi

[ -s ${outdir}/kmer_flanktooshort_4rm.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmer_flanktooshort_4rm.txt ${outdir}/preprocssing
fi

[ -s ${outdir}/kmer_flanktooshort_flkcoor.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmer_flanktooshort_flkcoor.txt ${outdir}/preprocssing
fi

[ -s ${outdir}/kmer_forblast.fasta ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmer_forblast.fasta ${outdir}/preprocssing
fi

[ -s ${outdir}/sigk_noN.fasta ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/sigk_noN.fasta ${outdir}/preprocssing
fi

[ -s ${outdir}/sigk_withN.fasta ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/sigk_withN.fasta ${outdir}/preprocssing
fi

mkdir ${outdir}/kmers_withN

[ -s ${outdir}/kmer_with_missinggenomes.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmer_with_missinggenomes.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/kmer_genomeappearonce.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmer_genomeappearonce.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/kmer_with_multi_hits.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmer_with_multi_hits.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/kmer_with_align_issue.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmer_with_align_issue.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/kmer_with_align_len.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmer_with_align_len.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/kmer_with_ID_E_issue_k.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmer_with_ID_E_issue_k.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/filterk_out_summary.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/filterk_out_summary.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/myout.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myout.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/myflk_behave_pheno.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myflk_behave_pheno.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/mysplitk_out.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/mysplitk_out.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/rows_for_process.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/rows_for_process.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/myshort_splitk_out_uniq.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myshort_splitk_out_uniq.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/myintactkwithN_out.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myintactkwithN_out.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/myintactkwithN_rev0fwd1_kmer4plot.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myintactkwithN_rev0fwd1_kmer4plot.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/myintactkwithN_rev0fwd1.png ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myintactkwithN_rev0fwd1.png ${outdir}/kmers_withN
fi

[ -s ${outdir}/myintactkwithN_rev0fwd1_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myintactkwithN_rev0fwd1_set.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/myintactkwithN_rev1fwd0_kmer4plot.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myintactkwithN_rev1fwd0_kmer4plot.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/myintactkwithN_rev1fwd0.png ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myintactkwithN_rev1fwd0.png ${outdir}/kmers_withN
fi

[ -s ${outdir}/myintactkwithN_rev1fwd0_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myintactkwithN_rev1fwd0_set.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/myintactkwithN_other_kmer4plot.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myintactkwithN_other_kmer4plot.txt ${outdir}/kmers_withN
fi

[ -s ${outdir}/myintactkwithN_other.png ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myintactkwithN_other.png ${outdir}/kmers_withN
fi

[ -s ${outdir}/myintactkwithN_other_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myintactkwithN_other_set.txt ${outdir}/kmers_withN
fi

mkdir ${outdir}/kmers_noN

[ -s ${outdir}/kmernoN_length.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmernoN_length.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/kmer_with_missinggenomes_noN.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmer_with_missinggenomes_noN.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/kmer_with_multi_hits_noN.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmer_with_multi_hits_noN.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/kmer_with_align_len_noN.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmer_with_align_len_noN.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/kmer_with_ID_E_issue_noN.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/kmer_with_ID_E_issue_noN.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/filterk_out_summary_noN.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/filterk_out_summary_noN.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/mynoN_out.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/mynoN_out.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/rows_for_process_NoN.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/rows_for_process_NoN.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/myflk_behave_pheno_NoN.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myflk_behave_pheno_NoN.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/myNoNintactk_out.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myNoNintactk_out.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/myNoNintactk_rev0fwd1_kmer4plot.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myNoNintactk_rev0fwd1_kmer4plot.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/myNoNintactk_rev0fwd1.png ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myNoNintactk_rev0fwd1.png ${outdir}/kmers_noN
fi

[ -s ${outdir}/myNoNintactk_rev0fwd1_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myNoNintactk_rev0fwd1_set.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/myNoNintactk_rev1fwd0_kmer4plot.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myNoNintactk_rev1fwd0_kmer4plot.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/myNoNintactk_rev1fwd0.png ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myNoNintactk_rev1fwd0.png ${outdir}/kmers_noN
fi

[ -s ${outdir}/myNoNintactk_rev1fwd0_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myNoNintactk_rev1fwd0_set.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/myNoNintactk_other_kmer4plot.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myNoNintactk_other_kmer4plot.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/myNoNintactk_other.png ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myNoNintactk_other.png ${outdir}/kmers_noN
fi

[ -s ${outdir}/myNoNintactk_other_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myNoNintactk_other_set.txt ${outdir}/kmers_noN
fi

[ -s ${outdir}/splitk_plots ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/splitk_plots ${outdir}/kmers_withN
fi

mkdir ${outdir}/allintack_combinedplots
[ -s ${outdir}/myallintactk_rev0fwd1.png ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myallintactk_rev0fwd1.png ${outdir}/allintack_combinedplots
fi

[ -s ${outdir}/myallintactk_rev1fwd0.png ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myallintactk_rev1fwd0.png ${outdir}/allintack_combinedplots
fi

[ -s ${outdir}/myallintactk_rev0fwd1_kmer4plot.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myallintactk_rev0fwd1_kmer4plot.txt ${outdir}/allintack_combinedplots
fi

[ -s ${outdir}/myallintactk_rev0fwd1_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myallintactk_rev0fwd1_set.txt ${outdir}/allintack_combinedplots
fi

[ -s ${outdir}/myallintactk_rev1fwd0_kmer4plot.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myallintactk_rev1fwd0_kmer4plot.txt ${outdir}/allintack_combinedplots
fi

[ -s ${outdir}/myallintactk_rev1fwd0_set.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
mv ${outdir}/myallintactk_rev1fwd0_set.txt ${outdir}/allintack_combinedplots
fi

[ -s ${outdir}/myNoNintactk_rev0fwd1_set_no1stline.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
rm ${outdir}/myNoNintactk_rev0fwd1_set_no1stline.txt
fi

[ -s ${outdir}/myNoNintactk_rev1fwd0_set_no1stline.txt ] && present=1 || present=0
if [[ ${present} -eq 1 ]]
then
rm ${outdir}/myNoNintactk_rev1fwd0_set_no1stline.txt
fi

gzip ${gen/.gz}
