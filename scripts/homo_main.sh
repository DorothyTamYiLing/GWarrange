#Usage: bash script/homo_main.sh -gff example_data/C505_NZ_CP011687.1.gff -fna example_data/C505_NZ_CP011687.1.fna 


#run this script inside first level directory, gff and genome should be in example_data
#set default parameters in getopts
thread_blast=8
freq=2
idcov="80_80"
dist=1000

while test $# -gt 0; do
           case "$1" in
                -thread_blast)
                    shift
                    thread_blast=$1
                    shift
                    ;;
                -gff)
                    shift
                    gff=$1
                    shift
                    ;;
                -fna)
                    shift
                    fna=$1
                    shift
                    ;;
                -freq)
                    shift
                    freq=$1
                    shift
                    ;;
                -idcov)
                    shift
                    idcov=$1
                    shift
                    ;;
                -dist)
                    shift
                    dist=$1
                    shift
                    ;;
                *)
                   echo "$1 is not a recognized flag!"
                   return 1;
                   ;;
          esac
  done  

echo "gff file : $gff";
echo "ref genome fasta file : $fna";
echo "homo freq : $freq";
echo "blast id_coverage: $idcov";
echo "number of threads in blast: $thread_blast";
echo "maximum distance to be defined as sam cluster: $dist";



#gff="example_data/PROKKA_sim1.gff"
#fna="example_data/sim1.fasta"

#remove top few lines containing "#"
grep -v "#" ${gff} > tmpfile && mv tmpfile ${gff}

#slice out each seqeunce in gff from the genome, output file: homo_gffseq_frgenome.fasta
python script/homo_cutgffseq_frgenome.py --gff ${gff} --genome ${fna}

#cut out blast id identity and query coverage
homo_id=$(echo ${idcov} | cut -d "_" -f1)
homo_cov=$(echo ${idcov} | cut -d "_" -f2)

#blast gff seq with genome
blastn -query homo_gffseq_frgenome.fasta -subject ${fna} -perc_identity ${homo_id} -qcov_hsp_perc ${homo_cov} -outfmt 6 -out homo_blastgffgen_out.txt
#blastn -query homo_gffseq_frgenome.fasta -subject ${fna} -outfmt 6 -out homo_blastgffgen_out.txt


#keep the gff seq with set blast hit freq , outout : homo_freq_list.txt
Rscript script/homo_freqliist.R --freq ${freq}

#extract seq in homo_freq_list.txt from homo_gffseq_frgenome.fasta
python3 script/getfastafrlist.py --list homo_freq_list.txt --input homo_gffseq_frgenome.fasta --out homo_freq.fna


#blast homo_freq.fna with itself
blastn -perc_identity ${homo_id} -qcov_hsp_perc ${homo_cov} -query homo_freq.fna -subject homo_freq.fna -outfmt 6 -out homo_homo_blastout.txt

#deduplicate homo_freq.fna based on blast output, output file : homo_deduplist.txt
#echo "deduplicate homo_freq"
Rscript script/homo_dedup.R

#extract seq in homo_deduplist.txt from gffseq_frgenome.fasta
python3 script/getfastafrlist.py --list homo_deduplist.txt --input homo_gffseq_frgenome.fasta --out homodedup.fna


#blast homodedup.fna with genome
blastn -query homodedup.fna -perc_identity ${homo_id} -qcov_hsp_perc ${homo_cov} -subject ${fna} -outfmt 6 -out homo_blastdedup_out.txt
#blastn -query homodedup.fna  -subject ${fna} -outfmt 6 -out homo_blastdedup_out.txt

#determine unique hoiomo cluters in genome
#echo "determine unique hoiomo cluters in genome"
Rscript script/homo_cluster_size.R --dist ${dist}

#make occurrence file
grep ">" homodedup.fna | sed 's/>//g;s/ ID=/delimiterID=/g;s/ /_/g;s/delimiter/\t/g' > header.txt
join -1 1 -2 1 -o 2.1,1.2,2.2 -e 0 homo_occurence.txt header.txt > homo_occurence_1.txt 
echo "homo_id    occurrence    annotation" | cat - homo_occurence_1.txt  > homo_occurence.txt

#remove intermediate files
rm homo_occurence_1.txt
rm homo_blastdedup_out.txt
rm homo_blastgffgen_out.txt
rm homo_freq_list.txt
rm homo_freq.fna
rm homo_gffseq_frgenome.fasta
rm homo_homo_blastout.txt

mkdir output_homo

mv homo* output_homo


