#set default parameters in getopts
extend=5
merge=3
switch="on"

while getopts g:o:a:i:e:m:s: flag
do
    case "${flag}" in
        g) gen=${OPTARG};;
        o) outdir=${OPTARG};;
        a) fixgene=${OPTARG};;
        i) myIS=${OPTARG};;
        e) extend=${OPTARG};;
        m) merge=${OPTARG};;
        s) switch=${OPTARG};;
    esac
done
echo "gen: $gen";
echo "outdir: $outdir";
echo "first gene in genome: $fixgene";
echo "IS fasta sequence: $myIS";
echo "IS extend: $extend";
echo "IS merge: $merge";
echo "switch for turning on/off minimal ISmerge: $switch";


#gunzip the genome set 
gunzip $gen

#create output directory
mkdir $outdir

#blast the first gene with genome 
#blastn -query ./genome_rearrangement/first_gene/gidA.fasta -subject test.fasta -outfmt 6 -perc_identity 80 -qcov_hsp_perc 80 -evalue 1e-5 -out blastfirstgene_out.txt
blastn -query ${fixgene} -subject ${gen/.gz} -outfmt 6 -perc_identity 80 -qcov_hsp_perc 80 -evalue 1e-5 -out ${outdir}/blastfirstgene_out.txt

#fix the genome to same orientation according to specified gene
python3 fix_genome.py --input ${gen/.gz} --coor ${outdir}/blastfirstgene_out.txt --outdir ${outdir}

#blast the IS481 with genome 
#blastn -query TOHAMA1_IS481_27283to28335.fasta -subject clus1clus2_47.fna -outfmt 6 -out blastIS_out.txt
blastn -query ${myIS} -subject ${outdir}/fixed_genomes.fasta -outfmt 6 -out ${outdir}/blastIS_out.txt

#merge and extending IS
Rscript merge_IS.R --input ${outdir}/blastIS_out.txt --extend $extend --merge $merge --outdir ${outdir}

rm -r ${outdir}/ext${extend}_merge${merge}_ISreplaced_genomes

python3 iSreplace_2col.py --input ${outdir}/fixed_genomes.fasta  --coor ${outdir}/ext${extend}_merge${merge}_mergedIS.txt --out ${outdir}/ext${extend}_merge${merge}_ISreplaced_genomes

if [[ $merge -gt 3 && $switch == "on" ]]  #if merge arguments is provided and switch for minimal ISmerg is on as default
then
echo "also perform merging overlapping IS only"

Rscript merge_IS.R --input ${outdir}/blastIS_out.txt --extend $extend --merge 3

rm -r ${outdir}/ext${extend}_merge3_ISreplaced_genomes

python3 iSreplace_2col.py --input ${outdir}/fixed_genomes.fasta  --coor ${outdir}/ext${extend}_merge3_mergedIS.txt --out ${outdir}/ext${extend}_merge3_ISreplaced_genomes
fi
