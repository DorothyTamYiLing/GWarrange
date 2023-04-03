#set default parameters in getopts
flk_len=30
flk_dist=200000

while getopts k:g:p:l:f:d:o: flag
do
    case "${flag}" in
        k) sigk=${OPTARG};;
        g) gen=${OPTARG};;
        p) pheno=${OPTARG};;
        l) k_len=${OPTARG};;
        f) flk_len=${OPTARG};;
        d) flk_dist=${OPTARG};;
        o) outdir=${OPTARG};;
    esac
done
echo "sigk: $sigk";
echo "gen: $gen";
echo "pheno: $pheno";
echo "k_len: $k_len";
echo "flk_len: $flk_len";
echo "flk_dist: $flk_dist";
echo "outdir: $outdir";

#create output directory
mkdir $outdir

#python3 class_k.py --input ./example_data/clus1clus2_sigk_withN.fasta --outdir ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir
python3 class_k.py --input $sigk --outdir $outdir

gunzip $gen

#processing sig kmers with N, only when the input fasta file is present
[ -s ${outdir}/sigk_withN.fasta ] && withN=1 || withN=0

if [[ ${withN} -eq 1 ]]
then
echo "there is sigk withN"
bash filtering_kmer_and_blast.sh ${outdir}/sigk_withN.fasta ${gen/.gz} $outdir $flk_len
echo "now run make_flank_summary.R"
Rscript make_flank_summary.R --k.len $k_len --pheno $pheno --outdir $outdir --flkdist $flk_dist
else
echo "no withN sig k"
fi

#processing sig kmers without N, only when the input fasta file is present
[ -s ${outdir}/sigk_noN.fasta ] && noN=1 || noN=0

if [[ ${noN} -eq 1 ]]
then
echo  "there is sigk_noN.fasta"
blastn -query ${outdir}/sigk_noN.fasta -subject ${gen/.gz} -outfmt 6 -out ${outdir}/mynoN_out.txt
echo "no run process_sigkNoN.R"
Rscript process_sigkNoN.R --pheno $pheno --outdir $outdir
else
echo "no NoN sig k"
fi

gzip ${gen/.gz}
