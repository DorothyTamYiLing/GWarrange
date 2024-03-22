Extra requirements for the tutorials:
pyseer 1.3.10, 
unitig-caller 1.3.0, 
Frequency-based String Mining (lite)

Tip 1: To avoid files confusion, before running a new tutorial, it is advisable to remove all output files/folders from previous tutorial runs.
Tip 2: gff file of selected reference genome for estimating homologous sequence clusters can be converted from genebank file (.gbk, downloaded from NCBI) using online tools such as http://genome2d.molgenrug.nl/g2d_tools_conversions.html

![diagrams_flowchart](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/6c05fe4f-5092-466e-b974-a1384dff2b80)

Genome rearrangement pipeline summary chart


# Tutorial 1

This tutorial is based a subset of _Bordetella pertussis_ genomes as described in Weigand _et al_. 2019, in which chromosome structures are defined by MAUVE exhaustive pairwise alignment. A subset of 47 genomes that display two different chromosome structures (18 genomes with structure "1" and 29 genomes with structure "0" (See figure below) are used in GWAS with structure information as phenotype. Structure phenotype of two pairs of genomes are swapped for demonstration purpose.  

<img width="652" alt="Screenshot 2024-01-10 at 18 17 39" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/2ccf91f0-17bc-4b01-b8a1-b3744eb2f749">

Fig : Two different chromosome structures that are found among 47 _Bordetella pertussis_ genomes. 

##

Go to the top level of /genome_rearrangement directory
```
cd /path/to/genome_rearrangement
```

First, the genome assemblies multifasta file is prepared by concatenating the genome fasta files.
```
#concatenating genome fasta files for use
cat example_data/example_genomes/clus1clus2_47_genomes/*fasta.gz > example_data/clus1clus2_47.fna.gz
```

Genomes assemblies from which genome rearrangements are detected are re-orientated by a chosen gene. In the case of _Bordetella pertussis_, the gene is _gidA_ since it is the first gene after origin of replication. The location and orientation of _gidA_ in the genomes are obtained by aligning it with multifasta file of genome assemblies by BLAST.
```
#unzip the genome file if neccesasry
gunzip example_data/clus1clus2_47.fna.gz

#blast gidA with genomes
blastn -query example_data/gidA.fasta \
-subject example_data/clus1clus2_47.fna \
-outfmt 6 -out clus1clus2_47_gidA_out.txt
```

Then, genome assemblies are re-orientated according to the position and orientation of _gidA_ in the genomes, using the script fix_genome.py:
```
python3 scripts/fix_genome.py --input example_data/clus1clus2_47.fna --mycoor clus1clus2_47_gidA_out.txt
```

The output multifasta file name for the genomes in the same orientation is "fixed_genomes.fasta".

Genome rearrangrements in _Bordetella pertussis_ are believed to be largely mediated by homologous recombination between insertion sequence elements (IS elements), such as IS481 and IS110. Location of target IS elements in the genomes are obtained through BLAST. Sequences of more than one target IS element can be placed in the same multifasta file for obtaining their genome locations in all genomes at once.
```
blastn -query example_data/IS_NZ_CP025371.1.fasta \
-subject fixed_genomes.fasta \
-outfmt 6 -out clus1clus2_47_blastIS_out.txt
```

In addition, genome rearrangements in _Bordetella pertussis_ have also been observed to be mediated by homologous recombination of sequence blocks that consist of one or more IS element. These duplicated sequence blocks are found throughout the genome and can be as large as several thousands base pairs in size. To ensure sensitivity in genome rearrangement detection, it is advised to replace these sequence blocks **completely** with shorter placeholder sequence. Without additional information of the actual size of the homologous sequence blocks, sequences extending several thousands base pairs to both directions from each IS element can be replaced.

Here, sequences extending 7000bp to both directions from each IS element are replaced. IS elements that are no more than 200bp apart (after extension) in each genome are also "merged". Then, each of these "extended and merged" IS elements are replaced with shorter placedholder sequences (N x 15). A seperate set of IS-replaced genomes are also produced by enabling performing minimal IS extension (i.e. 100bp to both directions) and merging overlapping IS elements only (i.e. IS elements that are less than 3bp apart) through passing a string argument "on" to the -s flag.
```
bash scripts/merge_replace_IS.sh -g fixed_genomes.fasta -i clus1clus2_47_blastIS_out.txt -e 7000 -m 200 -s "on"
```

Each set of IS-replaced genomes using different IS merging and extending parameters are output into new directories "ext7000_merge200_ISreplaced_genomes" and "ext100_merge3_ISreplaced_genomes" respectively.

Statistics on the sizes and distances between each pair of adjacent "extended and merged" IS are printed in files *_mergedISstat.txt.

Prior to GWAS, each set of IS-replaced genomes using different IS merging and extending parameters are used for generating kmers.

```
#For ext7000_merge200_ISreplaced_genomes set
cd ext7000_merge200_ISreplaced_genomes

#generating fsm-ite input file
for f in *_ext7000_merge200_ISreplaced.fasta; do id=$(basename "$f" _ext7000_merge200_ISreplaced.fasta); echo $id $f; done > clus1clus2_47_input.list

#generating kmers with size of 200 bases with minor allele frequency 0.05
fsm-lite -l clus1clus2_47_input.list -v -s 3 -S 44 -t tmp -m 200 -M 200 | gzip - > clus1clus2_47_ext7000merge200_k200_output.txt.gz

########################################################################

#For ext100_merge3_ISreplaced_genomes set
cd ext100_merge3_ISreplaced_genomes

for f in *_ext100_merge3_ISreplaced.fasta; do id=$(basename "$f" _ext100_merge3_ISreplaced.fasta); echo $id $f; done > clus1clus2_47_input.list

#generating kmers with size of 200 bases with minor allele frequency 0.05
fsm-lite -l clus1clus2_47_input.list -v -s 3 -S 44 -t tmp -m 200 -M 200 | gzip - > clus1clus2_47_ext100merge3_k200_output.txt.gz

```

Then, a kmer-based GWAS is conducted using pyseer with an aim to identify kmers whose presence-absence patterns are associated with chromosome structures phenotype. Population structure is not controlled.

```
#Run inside corresponding *_ISreplaced_genomes directory

#adding header to phenotype file for pyseer input format
echo "samples binary" | cat - ../example_data/clus1clus2_pheno.txt > ../example_data/clus1clus2_pheno_4pyseer.txt

#run pyseer

#For ext7000_merge200_ISreplaced_genomes set
pyseer --phenotypes ../example_data/clus1clus2_pheno_4pyseer.txt \
--kmers clus1clus2_47_ext7000merge200_k200_output.txt.gz \
--no-distances \
--min-af 0.05 --max-af 0.95 \
--print-samples --output-patterns kmer_patterns.txt \
> clus1clus2_47_ext7000merge200_k200_MAF0.05_nopopctrl

#For ext100_merge3_ISreplaced_genomes set
pyseer --phenotypes ../example_data/clus1clus2_pheno_4pyseer.txt \
--kmers clus1clus2_47_ext100merge3_k200_output.txt.gz \
--no-distances \
--min-af 0.05 --max-af 0.95 \
--print-samples --output-patterns kmer_patterns.txt \
> clus1clus2_47_ext100merge3_k200_MAF0.05_nopopctrl

```

Generate the number of unique patterns and p value significance threshold information:
```
#Run inside corresponding *_ISreplaced_genomes directory

python3 ../scripts/count_patterns.py kmer_patterns.txt > count_pattern.txt

#count_patterns.py is a script from pyseer package for calculating p-value threshold using Bonferroni correction

```
Extract kmers with p value below the significance threshold:
```
#Run inside corresponding *_ISreplaced_genomes directory

#For ext7000_merge200_ISreplaced_genomes set, threshold may vary slightly between different runs
awk '{ if ($4 <= 4.50E-04) { print } }' clus1clus2_47_ext7000merge200_k200_MAF0.05_nopopctrl > sigk_pyseer.txt

#For ext100_merge3_ISreplaced_genomes set, threshold may vary slightly between different runs
awk '{ if ($4 <= 5.43E-04) { print } }' clus1clus2_47_ext100merge3_k200_MAF0.05_nopopctrl > sigk_pyseer.txt
```

The sequences of kmers that are found to be significantly associated with structural phenotype are extracted and placed in a multifasta file.

Extract significant kmer sequences and convert them into fasta format:
```
#Run inside corresponding *_ISreplaced_genomes directory

#get the seqeunce only
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

paste -d \\n header.txt sigk_seq.txt > sigk_seq.fasta
```

Then, these kmers are blasted with the original genome set for studying potential genome rearrangements that are captured by them, implemented by the following script:

```
#run in the top level of /genome_rearrangement directory
cd /path/to/genome_rearrangement

#For ext7000_merge200_ISreplaced_genomes set
bash scripts/main.sh -k ext7000_merge200_ISreplaced_genomes/sigk_seq.fasta \
-g example_data/clus1clus2_47.fna.gz \
-p example_data/clus1clus2_pheno.txt -d 110000 -f 30 \
-o clus1clus2_47_ext7000_merge200_outdir -s 4300 -x 2 -y 1000

#For ext100_merge3_ISreplaced_genomes set
bash scripts/main.sh -k ext100_merge3_ISreplaced_genomes/sigk_seq.fasta \
-g example_data/clus1clus2_47.fna \
-p example_data/clus1clus2_pheno.txt -d 5000 -f 30 \
-o clus1clus2_47_ext100_merge3_outdir -s 4300 -x 2 -y 1000

```

Note that the value used for -d parameter should be larger than the "Maximum size of merged ISs" value in the corresponding *mergedISstat.txt file.

**Visualising genome rearrangements that are captured by kmer**

From genome set with 7000bp extension and 200bp merging, split kmers are found (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes. They can be found in clus1clus2_47_ext7000_merge200_outdir/kmers_withN/mysplitk_out.txt.

1) Plotting split kmers for visualising rearrangement boundaries

Since kmers contain highly redundant information, only kmers with unique information, i.e. unique behaviour count and proportion in case/control genomes, unique genome position (represented by the mean StartL when k-mers are intact in case genomes, rounded off to two significant digits, as indicated by the x flag), and unique forward/reverse intact k-mers count are kept. They can be found in output file clus1clus2_47_ext7000_merge200_outdir/kmers_withN/myshort_splitk_out_uniq.txt.

Four rearrangement boundaries are found, and they potentially refer to two inversion events, i.e. between 43000bp and 3600000bp, as well as between 1500000bp and 2500000bp, one inversion nested within the other. The four boundaries can be indicated by sixteen different significant split kmers that are mapped to each of the boundary, split in case/control genomes, and in forward/reverse orientation (plots of four split kmers are shown below as examples). Full information of these kmers can be found in output file clus1clus2_47_ext7000_merge200_outdir/kmers_withN/mysplitk_out.txt. Plots for split kmers can be found in folder clus1clus2_47_ext7000_merge200_outdir/kmers_withN/splitk_plots.

Inversion within genome region 43000 and 3600000, 43000bp boundary, kmer being intact in case genomes and split in control genomes, in forward orientation:

![kmer99_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/b5d132a7-8e65-453d-9b54-ffeb80334bc4)

Inversion within genome region 43000 and 3600000, 3600000bp boundary, kmer being intact in control genomes and split in case genomes, in reverse orientation:

![kmer933_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/bf62eeb5-ee75-449b-aab1-be13ddd4c77a)

Inversion within genome region 1500000 and 2500000, 1500000bp boundary, kmer being intact in control genomes and split in case genomes, forward orientation kmer:

![kmer932_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/0f97454f-107d-4643-881f-53f1fc1b7e45)

Inversion within genome region 1500000 and 2500000, 2500000bp boundary, kmer being intact in case genomes and split in control genomes, reverse orientation kmer:

![kmer931_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/531b0035-976e-4fca-bc53-ad21bbe5c6aa)


Height of arrows corresponds to proportion of case/control genomes.

2) Plotting intact kmers without N for visualising interior sequence content of rearrangements :

Genome position of intact kmers without N from /ext100_merge3_ISreplaced_genomes (minimal IS extension and merging overlapping IS only) are plotted. Only kmers with unqiue genome position information (by rounding off to the nearest multiple of 1000) are kept for plotting (as shown in clus1clus2_47_ext100_merge3_outdir/kmers_noN/*kmer4plot.txt files). Plots for intact kmers can be found in folder clus1clus2_47_ext100_merge3_outdir/kmers_noN.

Plot of intact kmers that show rearrangements in two genome regions that are significantly associated with structural phenotype.

![diagrams_17](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/fe65d012-960b-457a-818c-885bc71ded60)


(Left: intact kmers that are in forward orientation in majority of structure "1" genomes, as well as in reverse orientation in majority of structure "0" genomes;
Right: intact kmers that are in reverse orientation in majority of structure "1" genomes, as well as in forward orientation in majority of structure "0" genomes)

**Important notes:**

Inversion boundaries at 1500000bp and 2500000bp are not detected using ext100_merge3_ISreplaced_genomes. It is because homologous sequence blocks/IS elements clusters at these boundaries are not completely replaced by shorter placeholder sequences. This leads to the presence of homologous sequence in flanking sequences in the kmers, hence reducing sensitivity in detecting rearrangement boundaries. To ensure complete replacement of homologous sequence blocks/IS elements clusters, we can extend the genome coordinates of each repeated sequence for a number of base pairs in both directions, and/or to merge repeated sequences that are less than a number of base pairs distance apart.

Reference: Weigand, M.R., Williams, M.M., Peng, Y., Kania, D., Pawloski, L.C., Tondella, M.L. and CDC Pertussis Working Group, 2019. Genomic survey of Bordetella pertussis diversity, United States, 2000–2013. Emerging infectious diseases, 25(4), p.780.


# Tutorial 2

_Enterococcus faecium_'s genomes are known to be enriched with IS elements, which could play important role in their genome structure's diversification (Leavis _et al._ 2007). Genome structure of 75 _Enterococcus faecium_ were characterised by socru (Page _et al._ 2020), Among which, a subset of 32 genomes displaying two different chromosome structures (21 genomes with structure "0" and 11 genomes with structure "1") (See Fig. below) are used in GWAS, with structure information as phenotype. Structure phenotype of two pairs of genomes are swapped for demonstration purpose. 

<img width="680" alt="Screenshot 2024-01-10 at 18 17 45" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/a3791047-f56b-4565-85bc-43a7f1d897c2">

Fig : Two different chromosome structures found among 32 _Enterococcus faecium_ genomes. 

First, genomes asemblies from which genome rearrangements are detected are re-orientated by a chosen gene, i.e. dnaA. The location and orientation of dnaA in the genomens are obtained by aligning it with multifasta file of genome assemblies through BLAST.

Go to the top level of /genome_rearrangement directory
```
cd /path/to/genome_rearrangement
```
First, the genome assemblies multifasta file is prepared by concatenating genome fasta files.
```
#concatenating genome fasta files for use
cat example_data/example_genomes/Efaecium_32genomes/*fasta.gz > example_data/32genomes.fna.gz
```

Genomes assemblies from which genome rearrangements are detected are re-orientated by a chosen gene. In the case of _Enterococcus faecium_, the chosen gene can be dnaA. The location and orientation of dnaA in the genomes are obtained by blasting it with multifasta file of genome assemblies.

```
#unzip the genome file if neccesasry
gunzip example_data/32genomes.fna.gz

#blast dnaA with genomes
blastn -query example_data/dnaA.fasta \
-subject example_data/32genomes.fna \
-outfmt 6 -out 32genomes_dnaA_out.txt
```

Then, genome assemblies are re-orientated according to the position and orientation of dnaA in the genomes, using the script fix_genome.py:

```
python3 scripts/fix_genome.py --input example_data/32genomes.fna --mycoor 32genomes_dnaA_out.txt
```
The output file name for the genomes with same orientation is "fixed_genomes.fasta".

Location of IS elements in the genomes are obtained by aligning them with target IS elements sequences through BLAST. Sequences of more than one target IS elements can be placed in the same multifasta file for obtaining their genome locations in all genomes at once.

```
blastn -query example_data/IS30_IS1252_in_HOU503_657692to658645.fasta \
-subject fixed_genomes.fasta \
-outfmt 6 -out blastIS30_IS1252_in_HOU503_32genomes_out.txt
```

Here, sequences extending 7000bp to both directions from each IS elements are replaced. IS elements that are no more than 200bp apart (after extension) in each genome are also "merged". Then, each of these "extended and merged" IS elements are replaced with shorter placedholder sequences (N x 15). A seperate set of IS-replaced genomes are also produced by enabling performing minimal IS extension (i.e. 100bp) and merging overlapping IS only (i.e. IS that are less than 3 bp apart) through passing string argument "on" to the -s flag.
```
bash scripts/merge_replace_IS.sh -g fixed_genomes.fasta -i blastIS30_IS1252_in_HOU503_32genomes_out.txt -e 7000 -m 200 -s "on"
```
Prior to GWAS, each set of IS-replaced genomes using different merging and extending parameters are used for generating kmers and unitigs. Here, we use the kmer generation tool fsm-lite and unitig generation tool unitig-caller to generate kmers and unitigs from genomes in directory /ext100_merge3_ISreplaced_genomes.
```
cd ext100_merge3_ISreplaced_genomes

#generating fsm-ite input file
for f in *_ext100_merge3_ISreplaced.fasta; do id=$(basename "$f" _ext100_merge3_ISreplaced.fasta); echo $id $f; done > 32genomes_input.list

#generating kmers with size of 200 bases kmers present in at least 20 samples

fsm-lite -l 32genomes_input.list -v -s 20 -S 30 -t tmp -m 200 -M 200 | gzip - > ext100merge3_k200_output.txt.gz

#generating input file for unitig-caller
ls -d -1 $PWD/*.fasta > input.txt

#running unitig-caller
unitig-caller --call --rtab --pyseer --refs input.txt --out unitigcall_out

#fixing sample names in output file
sed -i 's/_ext100_merge3_ISreplaced//g' unitigcall_out.pyseer

gzip unitigcall_out.pyseer
```
Then, kmer-based and unitig-based GWAS are conducted using pyseer, with an aim to identify kmers/unitigs whose presence-absence patterns are associated with chromosome structure phenotype. Population structure is not controlled.
```
#adding header to phenotype file for pyseer input format
echo "samples binary" | cat - ../example_data/Efaecium32genomes_pheno_1swap.txt > ../example_data/Efaecium32genomes_pheno_1swap_4pyseer.txt

#run pyseer using k-mer and unitig file
pyseer --phenotypes ../example_data/Efaecium32genomes_pheno_1swap_4pyseer.txt \
--kmers ext100merge3_k200_output.txt.gz \
--no-distances \
--min-af 0.05 --max-af 0.95 \
--print-samples --output-patterns kmer_patterns.txt \
> ext100merge3_k200_min20samp_nopopctrl

pyseer --phenotypes ../example_data/Efaecium32genomes_pheno_1swap_4pyseer.txt \
--kmers unitigcall_out.pyseer.gz \
--no-distances \
--min-af 0.05 --max-af 0.95 \
--print-samples --output-patterns kmer_unitig_patterns.txt \
> ext100merge3_k200_min20samp_unitig_nopopctrl

```
Generate the number of unique patterns and p value significance threshold information:
```
python3 ../scripts/count_patterns.py kmer_patterns.txt > count_pattern.txt

python3 ../scripts/count_patterns.py kmer_unitig_patterns.txt > count_uni_pattern.txt
#count_patterns.py is a script from pyseer package for calculating p-value threshold using Bonferroni correction. 
```
Extract kmers/unitigs with p value below the significance threshold, threshold may vary slightly between different runs:
```
awk '{ if ($4 <= 2.29E-05) { print } }' ext100merge3_k200_min20samp_nopopctrl > sigk_pyseer.txt #kmer

awk '{ if ($4 <= 6.45E-06) { print } }' ext100merge3_k200_min20samp_unitig_nopopctrl > siguni_pyseer.txt #unitig
```
448,330 kmers and 4,853 unitigs are found to be significantly associated with chromosome structure. The sequences of which were extracted and placed in respective multifasta file.

Extract significant kmer sequences and convert them into fasta format
```
#get the seqeunce only
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

paste -d \\n header.txt sigk_seq.txt > sigk_seq.fasta
```
Due to the large number of significant kmers, only the kmers that contain "N" and the first 5000 kmers without "N" are used for structural analysis.
```
#classifying kmers into those containing "N" and those do not contain "N". Produce "sigk_noN.fasta" and "sigk_withN.fasta".
python3 ../scripts/class_k.py --input sigk_seq.fasta --outdir .

#extract first 5000 kmers without "N"
head -25000 sigk_noN.fasta > sigk_noN_5000.fasta

#combine the kmers with "N" and the first 5000 kmers without "N"
cat sigk_withN.fasta sigk_noN_5000.fasta > sigkwithN_noN5000.fasta
```

Extract significant unitig sequences and convert them into fasta format
```
#get the seqeunce only
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
	echo ">kmer""$c " >> header.txt
done

paste -d \\n header.txt siguni_seq.txt > siguni_seq.fasta
```

Then, these kmers and unitigs are blasted with the original genome set for studying potential genome rearrangement that are captured by them, implemented by the following script:


```
#run in the top level of /genome_rearrangement directory
cd /path/to/genome_rearrangement

#For ext100_merge3_ISreplaced_genomes set k-mers with N 
bash scripts/main.sh -k ext100_merge3_ISreplaced_genomes/sigkwithN_noN5000.fasta \
-g example_data/32genomes.fna.gz \
-p example_data/Efaecium32genomes_pheno_1swap.txt -d 3000 -f 30 \
-o Efaecium32genomes_ext100merge3_1swap_withNnoN5000_outdir -s 3000 -x 2 -y 1000

#For ext100_merge3_ISreplaced_genomes set unitigs
bash scripts/main.sh -k ext100_merge3_ISreplaced_genomes/siguni_seq.fasta \
-g example_data/32genomes.fna.gz \
-p example_data/Efaecium32genomes_pheno_1swap.txt -d 3000 -f 30 \
-o Efaecium32genomes_ext100merge3_1swap_unitig_outdir -s 3000 -x 2 -y 1000
```
Note that the value used for -d parameter should be larger than the "Maximum size of merged ISs" value in "ext100_merge3_mergedISstat.txt".

**Visualising genome rearrangements that are captured by kmer**

Split kmers are found (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes. They can be found in Efaecium32genomes_ext100merge3_1swap_withNnoN5000_outdir/kmers_withN/mysplitk_out.txt.

1) Plotting split kmers for visualising rearrangement boundaries

Since kmers contain highly redundant information, only kmers with unique information (genome position, case and control count and proportion) are kept. They can be found in output file Efaecium32genomes_ext100merge3_1swap_withNnoN5000_outdir/kmers_withN/myshort_splitk_out_uniq.txt.

Two rearrangement boundaries are found, and they potentially refer to a single inversion event between 72000bp and 2100000bp. All of the split kmers are found to be split in case geomes and intact in control genomes. This could be explained by the absence of IS elements in the rearrangement boundaries in case genomes, which has been confirmed by manual sequence check. The two boundaries can be indicated by four different significant split kmers that are mapped to each of the boundaries in forward/reverse orientation (plots of one of these kmers are shown below). Full information of these kmers can be found in output file mysplitk_out.txt. Plots for split kmers can be found in folder Efaecium32genomes_ext100merge3_1swap_withNnoN5000_outdir/kmers_withN/splitk_plots.

72000bp boundary, intact kmer in forward orientation:

![kmer98323_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/b012559f-bd0c-4ceb-b77e-53c0e07eadb5)

72000bp boundary, intact kmer in reverse orientation:
![kmer97292_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/237de0e6-6c80-4c5c-ae21-b9b68951af34)

2100000bp boundary, intact kmer in forward orientation:

![kmer99909_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/e27620d6-2ec5-4bd7-baaa-53f834bbae75)

2100000bp boundary, intact kmer in reverse orientation:

![kmer98897_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/4d33ac75-246d-473d-b189-45dc471db2c0)


No significant split unitigs are identified as unitig-caller are not able to generate unitigs containing placeholder sequences.

2) Plotting intact unitigs (without N) for visualising sequence content of rearrangement :

Genome positions of unitigs from /ext100_merge3_ISreplaced_genomes (minimal IS extension and merging overlapping IS only) are plotted. Only unitigs with unqiue genome position information (by rounding off to the nearest multiple of 1000) are kept for plotting (as shown in Efaecium32genomes_ext100merge3_1swap_unitig_outdir/kmers_noN/*kmer4plot.txt files). Plots for unitigs can be found in folder Efaecium32genomes_ext100merge3_1swap_unitig_outdir/kmers_noN.

Plot of unitigs that show rearrangements significantly associated with structure phenotype.

![diagrams_18](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/1d5f7076-33c3-4ccb-86f8-9ff4ee8e3cdd)


(Left: intact kmers that are in forward orientation in majority of structure "1" genomes, as well as in reverse orientation in majority of structure "0" genomes;
Right: intact kmers that are in reverse orientation in majority of structure "1" genomes, as well as in forward orientation in majority of structure "0" genomes)


Page, A.J., Ainsworth, E.V. and Langridge, G.C., 2020. socru: typing of genome-level order and orientation around ribosomal operons in bacteria. Microbial Genomics, 6(7).

Leavis, H.L., Willems, R.J.L., van Wamel, W.J.B., Schuren, F.H., Caspers, M.P.M. and Bonten, M.J.M., 2007. Insertion sequence–driven diversification creates a globally dispersed emerging multiresistant subspecies of E. faecium. PLoS pathogens, 3(1), p.e7.

# Tutorial 3

This tutorial is based on 468 _Bordetella pertussis_ genomes with pertactin (PRN) expression information. Among them, 165 genomes show the presence of pertactin expression and 303 show absence. Pertactin expression information is taken from the supplementary material summarised in Lefrancq _et al._ 2022.

Go to the top level of /genome_rearrangement directory
```
cd /path/to/genome_rearrangement
```

First, the genome assemblies multifasta file is prepared by concatenating genome fasta files.
```
#concatenating genome fasta files for use
cat example_data/example_genomes/PRN_468/*fasta.gz > example_data/PRN_468.fna.gz
```

Genomes assemblies from which genome rearrangements are detected are re-orientated by a chosen gene, i.e. gidA. The location and orientation of gidA in the genomens are obtained by aligning it with multifasta file of genome assemblies through BLAST.
```
#unzip the genome file if neccesasry
gunzip example_data/PRN_468.fna.gz

#blast gidA with genomes
blastn -query example_data/gidA.fasta \
-subject example_data/PRN_468.fna \
-outfmt 6 -out PRN_468_gidA_out.txt
```
Then, genome assemblies are re-orientated according to the position and orientation of gidA in the genomes, using the script fix_genome.py:
```
python3 scripts/fix_genome.py --input example_data/PRN_468.fna --mycoor PRN_468_gidA_out.txt
```
The output file name for the genomes in same orientation is "fixed_genomes.fasta".

Locations of IS elements in the genomes are obtained by aligning them with target IS sequences through BLAST. Sequences of more than one target IS element can be placed in the same multifasta file for obtaining their genome locations in all genomes at once.
```
blastn -query example_data/IS_NZ_CP025371.1.fasta \
-subject fixed_genomes.fasta \
-outfmt 6 -out PRN_468_blastIS_out.txt
```

Here, sequences extending 7000bp to both directions from each IS element are replaced. IS elements that are no more than 200bp apart (after extension) in each genome are also "merged". Then, each of these "extended and merged" IS element are replaced with shorter placedholder sequences (N x 15). A seperate set of IS-replaced genomes are also produced by enabling performing minimal IS extension (i.e. 100bp) and merging overlapping IS only (i.e. IS that are less than 3 bp apart) through passing string argument "on" to the -s flag.
```
bash scripts/merge_replace_IS.sh -g fixed_genomes.fasta -i PRN_468_blastIS_out.txt -e 7000 -m 200 -s "on"
```
Prior to GWAS, each set of IS-replaced genomes using different IS merging and extending parameters are used for generating kmers.
```
#For ext7000_merge200_ISreplaced_genomes set
cd ext7000_merge200_ISreplaced_genomes

#generating fsm-ite input file
for f in *_ext7000_merge200_ISreplaced.fasta; do id=$(basename "$f" _ext7000_merge200_ISreplaced.fasta); echo $id $f; done > PRN_468_input.list

#generating kmers with size of 200 bases with minor allele frequency 0.05

fsm-lite -l PRN_468_input.list -v -s 24 -S 444 -t tmp -m 200 -M 200 | gzip - > PRN_468_ext7000merge200_k200_output.txt.gz

##############################################################################

#For ext100_merge3_ISreplaced_genomes set
cd ext100_merge3_ISreplaced_genomes

#generating fsm-ite input file
for f in *_ext100_merge3_ISreplaced.fasta; do id=$(basename "$f" _ext100_merge3_ISreplaced.fasta); echo $id $f; done > PRN_468_input.list

#generating kmers with size of 200 bases with minor allele frequency 0.05
fsm-lite -l PRN_468_input.list -v -s 24 -S 444 -t tmp -m 200 -M 200 | gzip - > PRN_468_ext100merge3_k200_output.txt.gz

```
Then, a kmer-based GWAS is conducted using pyseer with an aim to identify kmers whose presence-absence patterns are associated with PRN expression phenotype. Population structure is controlled by phylogenetic similarity matrix.

```
#Run inside corresponding *_ISreplaced_genomes directory

#adding header to phenotype file for pyseer input format
echo "samples binary" | cat - ../example_data/prn_status_pheno.txt > ../example_data/prn_status_pheno_4pyseer.txt

#run pyseer

#For ext7000_merge200_ISreplaced_genomes set
pyseer --lmm --phenotypes ../example_data/prn_status_pheno_4pyseer.txt \
--kmers PRN_468_ext7000merge200_k200_output.txt.gz \
--similarity ../example_data/ClfML_kappa4.964_phylogeny_similarity.tsv \
--min-af 0.05 --max-af 0.95 --covariates ../example_data/covariates.txt --use-covariates 2 \
--print-samples --output-patterns kmer_patterns_covariate.txt \
> PRN468_ext7000merge200_k200_MAF0.05_covariate

#For ext100_merge3_ISreplaced_genomes set
pyseer --lmm --phenotypes ../example_data/prn_status_pheno_4pyseer.txt \
--kmers PRN_468_ext100merge3_k200_output.txt.gz \
--similarity ../example_data/ClfML_kappa4.964_phylogeny_similarity.tsv \
--min-af 0.05 --max-af 0.95 --covariates ../example_data/covariates.txt --use-covariates 2 \
--print-samples --output-patterns kmer_patterns_covariate.txt \
> PRN468_ext100merge3_k200_MAF0.05_covariate
```

Generate the number of unique patterns and p value significance threshold information (For both ext7000_merge200_ISreplaced_genomes set and ext100_merge3_ISreplaced_genomes set):
```
#Run inside corresponding *_ISreplaced_genomes directory

python3 ../scripts/count_patterns.py kmer_patterns_covariate.txt > count_pattern.txt

#count_patterns.py is a script from pyseer package for calculating p-value threshold using Bonferroni correction. 

```
Extract kmers with p value below the significance threshold (For both ext7000_merge200_ISreplaced_genomes set and ext100_merge3_ISreplaced_genomes set):
```
#Run inside corresponding *_ISreplaced_genomes directory

#For ext7000_merge200_ISreplaced_genomes set, threshold may vary slightly between different runs
awk '{ if ($4 <= 1.30E-04) { print } }' PRN468_ext7000merge200_k200_MAF0.05_covariate > sigk_pyseer.txt

#For ext100_merge3_ISreplaced_genomes set, threshold may vary slightly between different runs
awk '{ if ($4 <=  1.32E-04) { print } }' PRN468_ext100merge3_k200_MAF0.05_covariate > sigk_pyseer.txt
```

The sequences of kmers that are found to be significantly associated with the structure phenotype are extracted and placed in a multifasta file.

Extract significant kmer sequences and convert them into multifasta format (For both ext7000_merge200_ISreplaced_genomes set and ext100_merge3_ISreplaced_genomes set):
```
#Run inside corresponding *_ISreplaced_genomes directory

#get the seqeunce only
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

paste -d \\n header.txt sigk_seq.txt > sigk_seq.fasta
```

Then, these kmers are aligned through BLAST with the original genome set for studying potential genome rearrangements that are captured by them, implemented by the following script:

```
#run in the top level of /genome_rearrangement directory
cd /path/to/genome_rearrangement

#For ext7000_merge200_ISreplaced_genomes set
bash scripts/main.sh -k ext7000_merge200_ISreplaced_genomes/sigk_seq.fasta \
-g example_data/PRN_468.fna \
-p example_data/prn_status_pheno.txt -d 110000 -f 30 \
-o PRN_468_ext7000_merge200_outdir -s 4300 -x 2 -y 1000

#For ext100_merge3_ISreplaced_genomes set
bash scripts/main.sh -k ext100_merge3_ISreplaced_genomes/sigk_seq.fasta \
-g example_data/PRN_468.fna.gz \
-p example_data/prn_status_pheno.txt -d 5000 -f 30 \
-o PRN_468_ext100_merge3_outdir -s 4300 -x 2 -y 1000

```

Note that the value used for -d parameter should be larger than the "Maximum size of merged ISs" value in the corresponding *_mergedISstat.txt file.

**Visualising genome rearrangements that are captured by kmer**

1) Plotting split kmers for visualising rearrangement boundaries

No split kmer that indicated phenotype-associated rearrangement boundary is detected.

2) Plotting intact kmers without N for visualising interior sequence content of rearrangement :

Genome position of intact kmers without N from /ext100_merge3_ISreplaced_genomes (merging overlapping IS only) are plotted. Only kmers with unqiue genome position information (by rounding off to the nearest multiple of 1000) are used in the plot. Plots for intact kmers can be found in folder PRN_468_ext100_merge3_outdir/kmers_noN.

Plots of intact kmers that show interior rearranged sequence content that are significantly associated with PRN expression phenotype.

![diagrams_19](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/b94204fa-5bf6-41b1-85b0-fb89515e314b)

(Left: intact kmers that are in forward orientation in half of PRN+ (1) genomes, as well as reverse in majority of PRN- (0) genomes;
Right: intact kmers that are in reverse orientation in half of PRN+ (1) genomes, as well as forward in majority of PRN- (0) genomes)



**Important notes:**

Some of the significant intact kmers without placeholder sequence contain sequence of pertactin autotransporter (indicated by black arrows). These kmers are not found using ext7000_merge200_ISreplaced_genomes. This is because the gene pertactin autotransporter is located immediately next to an IS element in _pertusis_ genomes, and any genome rearrangement that sits completely within the "replaced IS" region will not be detected. 

Ref: Lefrancq, N., Bouchez, V., Fernandes, N., Barkoff, A.M., Bosch, T., Dalby, T., Åkerlund, T., Darenberg, J., Fabianova, K., Vestrheim, D.F. and Fry, N.K., 2022. Global spatial dynamics and vaccine-induced fitness changes of Bordetella pertussis. Science Translational Medicine, 14(642), p.eabn3253.

