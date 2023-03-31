#example of main.sh
#Usage: bash main.sh allsig_kmer_withN.fasta allsig_kmer_NoN.fasta 111_yearGWAS_genlist.fasta.gz  \
#/home/ubuntu/Dorothy/genome_rearrangement/phenotypes.tsv \
#/home/ubuntu/Dorothy/genome_rearrangement/output 200 30 2500

###define all the variables in the command###
#$1=k_input=allsig_kmer_withN.fasta  #with path
#$2=k_input=allsig_kmer_noN.fasta  #with path
#$3=gen_input=111_yearGWAS_genlist.fasta #with path
#$4=pheno=phenotypes.tsv  #with path
#$5=outdir=/home/ubuntu/Dorothy/genome_rearrangement/output
#$6=k_len=200
#$7=flnk_len=30
#$8=flkdist=2500
#############################################

ARG2=${2:-"no_NoN_intactk"}
ARG6=${6:-200}
ARG7=${7:-30}
ARG8=${8:-70000}

#create output directory
mkdir $5

#processing sig kmers with N, only when the input fasta file is present
if [[ $1!="no_withN_sigk" ]]
then
bash filtering_kmer_and_blast.sh $1 $3 $5 $ARG7

Rscript make_flank_summary.R --k.len $ARG6 --pheno $4 --outdir $5 --flkdist $ARG8
fi

#processing sig kmers without N, only when the input fasta file is present
if [[ $ARG2!="no_NoN_intactk" ]]
then
#blasting intact kmers with genomes 
blastn -query $2 -subject $3 -outfmt 6 -out $5/mynoN_out.txt

Rscript process_sigkNoN.R --pheno $4 --outdir $5
fi
