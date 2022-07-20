#Usage: bash main.sh allsig_kmer_withN.fasta 111_yearGWAS_genlist.fasta.gz  \
#/home/ubuntu/Dorothy/genome_rearrangement/phenotypes.tsv \
#/home/ubuntu/Dorothy/genome_rearrangement/output 200 30 2500


###define all the variables in the command
#$1=k_input=allsig_kmer_withN.fasta  #with path
#$2=gen_input=111_yearGWAS_genlist.fasta #with path
#$3=pheno=phenotypes.tsv  #with path
#$4=outdir=/home/ubuntu/Dorothy/genome_rearrangement/output
#$5=k_len=200
#$6=flnk_len=30
#$7=flkdist=2500

#variables in the pipeline (generated from input files)
k_num=$(grep ">" $1 | wc -l)
#k_num=$(grep ">" allsig_kmer_withN.fasta | wc -l)
gen_num=$(zcat $2 | grep ">" | wc -l)
#gen_num=$(zcat 111_yearGWAS_genlist.fasta.gz | grep ">" | wc -l)

#create output directory

#mkdir /home/ubuntu/Dorothy/genome_rearrangement/output
mkdir $4

bash filtering_kmer_and_blast.sh $1 $2 $4 $6
#bash filtering_kmer_and_blast.sh allsig_kmer_withN.fasta 111_yearGWAS_genlist.fasta.gz \
#/home/ubuntu/Dorothy/genome_rearrangement/output 30

#Rscript make_flank_summary.R --k.len 200 --pheno /home/ubuntu/Dorothy/genome_rearrangement/phenotypes.tsv --outdir /home/ubuntu/Dorothy/genome_rearrangement/output --flkdist 2500
Rscript make_flank_summary.R --k.len $5 --pheno $3 --outdir $4 --flkdist $7
