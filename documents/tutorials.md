
<img width="591" alt="Screenshot 2023-09-18 215805" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/f3501563-7c65-4c04-98c5-12a13ed3c3ad">

Genome rearrangement pipeline summary chart

Tip: To avoid files confusion, before start running a new tutorial, it is advisable to remove all output files/folders from previous tutorial run.

# Tutorial 1

This tutorial is based a subset of _Bordetella pertussis_ genomes as described in Weigand _et al_. 2019, in which chromosome structures are defined by MAUVE exhaustive pairwise alignment. A subset of 47 genomes that display two different chromosome structures (18 genomes with structure 1 and 29 genomes with structure 2) (See figure below) aree used in GWAS with structure information as phenotype. Structure phenotype of two pairs of genomes are swapped for demonstration purpose.  

<img width="541" alt="Screenshot 2023-10-05 234547" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/38a5f7e3-611b-4748-8a27-ab4c4e96af6b">

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

Genomes assemblies from which genome rearrangements are detected are re-orientated by a chosen gene. In the case of _Boredetella pertussis_, the gene is _gidA_ since it is the first gene after origin of replication. The location and orientation of _gidA_ in the genomes are obtained by blasting it with multifasta file of genome assemblies.
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

The output multifasta file name for the genomes with the same orientation is "fixed_genomes.fasta".

Genome rearrangments in _Bordetella pertussis_ are believed to be largely mediated by homologous recombination between insertion sequence (IS) elements (such as IS481 and IS110). Location of target IS elements in the genomes are obtained through BLAST. Sequences of more than one target IS element can be placed in the same multifasta file for obtaining their genome locations in all genomes at once.
```
blastn -query example_data/IS_NZ_CP025371.1.fasta \
-subject fixed_genomes.fasta \
-outfmt 6 -out clus1clus2_47_blastIS_out.txt
```

In addition, genome rearrangements in _Bordetella pertussis_ have also been observed to be mediated by homologous recombination of sequence blocks that consist of one or more IS elements. These duplicated sequence blocks are found throughout the genome and can be as large as several thousand base pairs in size. To ensure sensitivity in genome rearrangement detection, it is advised to replace these sequence blocks **completely** with shorter placeholder sequence. Without additional information of the actual size of the homologous sequence blocks, sequences extending several thousands base pairs to both direction from each IS can be replaced.

Here, sequences extending 7000bp to both direction from each IS are replaced. IS elements that are no more than 200bp apart (after extension) in each genome are also "merged". Then, each of these "extended and merged" IS elements are replaced with shorter placedholder sequences (N x 15). A seperate set of IS-replaced genomes are also produced by enabling performing minimal IS extension (i.e. 100bp to both direction) and merging overlapping IS only (i.e. ISs that are less than 3bp apart) through passing a string argument "on" to the -s flag.
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

Then, a kmer-based GWAS was conducted using pyseer with an aim to identify kmers whose presence-absence patterns were associated with chromosome structures phenotype. Population structure is not controlled.

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

Generate number of unique patterns and p value significance threshold information:
```
#Run inside corresponding *_ISreplaced_genomes directory

python3 ../scripts/count_patterns.py kmer_patterns.txt > count_pattern.txt

#count_patterns.py is a script from pyseer package for calculating p-value threshold using Bonferroni correction

```
Extract kmers with p value below the significance threshold:
```
#Run inside corresponding *_ISreplaced_genomes directory

#For ext7000_merge200_ISreplaced_genomes set
awk '{ if ($4 <= 4.59E-04) { print } }' clus1clus2_47_ext7000merge200_k200_MAF0.05_nopopctrl > sigk_pyseer.txt

#For ext100_merge3_ISreplaced_genomes set
awk '{ if ($4 <= 5.62E-04) { print } }' clus1clus2_47_ext100merge3_k200_MAF0.05_nopopctrl > sigk_pyseer.txt
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

From genome set with 7000bp extension and 200bp merging, 1008 kmers are found to be split (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes. They can be found in clus1clus2_47_ext7000_merge200_outdir/kmers_withN/mysplitk_out.txt.

1) Plotting split kmers for visualising rearrangement boundaries

Since kmers contain highly redundant information, only kmers with unique information (genome position, case and control count and proportion) are kept. They can be found in output file clus1clus2_47_ext7000_merge200_outdir/kmers_withN/myshort_splitk_out_uniq.txt.

Four rearrangement boundaries are found, and they potentially refer to two inversion events, i.e. between 43000bp and 3600000bp, as well as between 1500000bp and 2500000bp, one inversion nested within the other. The four boundaries can be indicated by sixteen different significant split kmers that are mapped to each of the boundary, split in case/control genomes, and in forward/reverse orientation (plots of four split kmers are shown below as examples). Full information of these kmers can be found in output file clus1clus2_47_ext7000_merge200_outdir/kmers_withN/mysplitk_out.txt. Plots for split kmers can be found in folder clus1clus2_47_ext7000_merge200_outdir/kmers_withN/splitk_plots.

Inversion within genome region 43000 and 3600000, 43000bp boundary, kmer being intact in case genomes and split in control genomes, in forward orientation:

<img width="707" alt="kmer940" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/4a747e7c-93e3-44b6-abe2-81b90e758370">

Inversion within genome region 43000 and 3600000, 3600000bp boundary, kmer being intact in control genomes and split in case genomes, in reverse orientation:

<img width="728" alt="kmer985" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/17d5a869-2cd7-49a9-884c-684d2a30ac17">

Inversion within genome region 1500000 and 2500000, 1500000bp boundary, kmer being intact in control genomes and split in case genomes, forward orientation kmer:

<img width="733" alt="kmer882" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/80de6ab4-952d-4c85-b15c-4cc74d80f990">

Inversion within genome region 1500000 and 2500000, 2500000bp boundary, kmer being intact in case genomes and split in control genomes, reverse orientation kmer:

<img width="740" alt="kmer881" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/4197c388-32b7-4456-acba-47a78f3ffdb8">

Height of arrows corresponds to proportion of case/control genomes.

2) Plotting intact kmers without N for visualising sequence content of rearrangement :

Genome position of intact kmers without N from /ext100_merge3_ISreplaced_genomes (minimal IS extension and merging overlapping IS only) are plotted. Only kmers with unqiue genome position information (by rounding off to the nearest multiple of 1000) are kept for plotting (as shown in clus1clus2_47_ext100_merge3_outdir/kmers_noN/*kmer4plot.txt files). Plots for intact kmers can be found in folder clus1clus2_47_ext100_merge3_outdir/kmers_noN.

Plot of intact kmers that show rearrangements in two genome regions that are significantly associated with structural phenotype.

![myNoNintactk_rev0fwd1](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/afde6c42-add7-4e07-ada9-7eedd348f3d3)

![myNoNintactk_rev1fwd0](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/b8a70fdd-1a4a-4cc0-b60e-e873a00a8d42)


(Above: 10 intact kmers that are in forward orientation in majority of structure "1" genomes, as well as in reverse orientation in majority of structure "0" genomes;
Below: 13 intact kmers that are in reverse orientation in majority of structure "1" genomes, as well as in forward orientation in majority of structure "0" genomes)

**Important notes:**

Inversion boundaries at 1500000bp and 2500000bp are not detected using ext100_merge3_ISreplaced_genomes. This is because the whole homologous sequence block/IS clusters at these boundaries are not completely replaced by shorter placeholder sequences, this leads to the presence of homologous sequence within flanking sequences in the kmers, hence reducing sensitivity in detecting rearrangement boundaries. To ensure replacing the whole homologous sequence block/IS clusters by short placeholder sequences, we can extend the genome coordinates of each repeated sequence for a number of base pairs in both directions, and/or to merge repeated sequences that are less than a number of base pairs distance apart.

Reference: Weigand, M.R., Williams, M.M., Peng, Y., Kania, D., Pawloski, L.C., Tondella, M.L. and CDC Pertussis Working Group, 2019. Genomic survey of Bordetella pertussis diversity, United States, 2000–2013. Emerging infectious diseases, 25(4), p.780.


# Tutorial 2

_Enterococcus faecium_'s genomes are known to be enriched with IS elements, which could play important role in their genome structure's diversification (Leavis _et al._ 2007). Genome structure of 75 _Enterococcus faecium_ were characterised by socru (Page _et al._ 2020), Among which, a subset of 32 genomes displaying two different chromosome structures (21 genomes with structure "0" and 11 genomes with structure "1") (See Fig. below) are used in GWAS, with structure information as phenotype. Structure phenotype of two pairs of genomes sre swapped for demonstration purpose. 

<img width="602" alt="Screenshot 2023-10-05 234635" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/05b4aa5b-71ea-4903-963f-ba5542adcac5">

Fig : Two different chromosome structures are found among 32 _Enterococcus faecium_ genomes. 

First, genomes asemblies from which genome rearrangements are detected are re-orientated by a chosen gene, i.e. dnaA. The location and orientation of dnaA in the genomens are obtained by blasting it with multifasta file of genome assemblies.

Go to the top level of /genome_rearrangement directory
```
cd /path/to/genome_rearrangement
```
First, the genome assemblies multifasta file is prepared by concatenating genome fasta files.
```
#concatenating genome fasta files for use
cat example_data/example_genomes/Efaecium_32genomes/*fasta.gz > example_data/32genomes.fna.gz
```

Genomes asemblies from which genome rearrangements are detected are re-orientated by a chosen gene. In the case of _Enterococcus faecium_, the chosen gene can be dnaA. The location and orientation of dnaA in the genomes are obtained by blasting it with multifasta file of genome assemblies.

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

Location of IS elements in the genomes are obtained by blasting with target IS sequences. Sequences of more than one target IS elements can be placed in the same multifasta file for obtaining their genome locations in all genomes at once.

```
blastn -query example_data/IS30_IS1252_in_HOU503_657692to658645.fasta \
-subject fixed_genomes.fasta \
-outfmt 6 -out blastIS30_IS1252_in_HOU503_32genomes_out.txt
```

Here, sequences extending 7000bp to both directions from each IS are replaced. IS elements that are no more than 200bp apart (after extension) in each genome are also "merged". Then, each of these "extended and merged" IS element are replaced with shorter placedholder sequences (N x 15). A seperate set of IS-replaced genomes are also produced by enabling performing minimal IS extension (i.e. 100bp) and merging overlapping IS only (i.e. IS that are less than 3 bp apart) through passing string argument "on" to the -s flag.
```
bash scripts/merge_replace_IS.sh -g fixed_genomes.fasta -i blastIS30_IS1252_in_HOU503_32genomes_out.txt -e 7000 -m 200 -s "on"
```
Prior to GWAS, each set of IS-replaced genomes using different IS merging and extending parameters are used for generating kmers and unitigs. Here, we use the kmer generation tool fsm-lite and unitig generation tool unitig-caller to generate kmers and unitigs from genomes in directory /ext100_merge3_ISreplaced_genomes.
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
```
Then, kmer-based and unitig-based GWAS are conducted using pyseer with an aim to identify kmers/unitigs whose presence-absence patterns are associated with chromosome structure phenotype. Population structure is not controlled.
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
Generate number of unique patterns and p value significance threshold information:
```
python3 ../scripts/count_patterns.py kmer_patterns.txt > count_pattern.txt

../scripts/count_patterns.py unitigs_kmer_patterns.txt > count_uni_pattern.txt
#count_patterns.py is a script from pyseer package for calculating p-value threshold using Bonferroni correction. 
```
Extract kmers/unitigs with p value below the the significance threshold:
```
awk '{ if ($4 <= 2.29E-05) { print } }' ext100merge3_k200_min20samp_nopopctrl > sigk_pyseer.txt #kmer

awk '{ if ($4 <= 7.62E-06) { print } }' ext100merge3_k200_min20samp_unitigs_nopopctrl > siguni_pyseer.txt #unitig
```
448,330 kmers and 3,737 unitigs are found to be significantly associated with chromosome structure. The sequences of which were extracted and placed in respective multifasta file.

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
Due to the large number of significant kmers, only the kmers that contain "N" and the first 5000 kmers without "N" were used for structural analysis.
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

rm header.txt   #remove any existing header file

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
-g example_data/Efaecium_32genomes.fna.gz \
-p example_data/Efaecium32genomes_pheno_1swap.txt -d 3000 -f 30 \
-o Efaecium32genomes_ext100merge3_1swap_withNnoN5000_outdir -s 3000 -x 2 -y 1000

#For ext100_merge3_ISreplaced_genomes set unitigs
bash scripts/main.sh -k ext100_merge3_ISreplaced_genomes/siguni_seq.fasta \
-g example_data/Efaecium_32genomes.fna.gz \
-p example_data/Efaecium32genomes_pheno_1swap.txt -d 3000 -f 30 \
-o Efaecium32genomes_ext100merge3_1swap_unitig_outdir -s 3000 -x 2 -y 1000
```
Note that the value used for -d parameter should be larger than the "Maximum size of merged ISs" value in "ext100_merge3_mergedISstat.txt".

**Visualising genome rearrangements that are captured by kmer**

215 kmers are found to be split (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes. They can be found in Efaecium32genomes_ext100merge3_1swap_withNnoN5000_outdir/kmers_withN/mysplitk_out.txt.

1) Plotting split kmers for visualising rearrangement boundaries

Since kmers contain highly redundant information, only kmers with unique information (genome position, case and control count and proportion) are kept. They can be found in output file Efaecium32genomes_ext100merge3_1swap_withNnoN5000_outdir/kmers_withN/myshort_splitk_out_uniq.txt.

Two rearrangement boundaries are found, and they potentially refer to a single inversion event between 72000bp and 2100000bp. All of the split kmers are found to be split in case geomes and intact in control genomes. This could be explained by the absence of IS elements in the rearrangement boundaries in case genomes, which has been confirmed by manual sequence check. The two boundaries can be indicated by four different significant split kmers that are mapped to each of the boundaries in forward/reverse orientation (plots of one of these kmers are shown below). Full information of these kmers can be found in output file mysplitk_out.txt. Plots for split kmers can be found in folder Efaecium32genomes_ext100merge3_1swap_withNnoN5000_outdir/kmers_withN/splitk_plots.

72000bp boundary, intact kmer in forward orientation:

<img width="659" alt="Screenshot 2023-10-17 153038" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/b67776fb-34ea-48a6-a118-3b4f6e830b51">

72000bp boundary, intact kmer in reverse orientation:

<img width="740" alt="Screenshot 2023-10-17 153257" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/6a15eedf-f51f-4343-b16f-40a1042d66a6">

2100000bp boundary, intact kmer in forward orientation:

<img width="814" alt="Screenshot 2023-10-17 153223" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/2376419c-8dc7-4257-b98e-3d3ca3717c82">


2100000bp boundary, intact kmer in reverse orientation:

<img width="710" alt="Screenshot 2023-10-17 153153" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/9de2d524-3c55-4550-9b16-e2bd529e82aa">

No significant split unitigs are identified as unitig-caller are not able to generate unitigs containing placeholder sequences.

2) Plotting intact unitigs (without N) for visualising sequence content of rearrangement :

Genome position of unitigs from /ext100_merge3_ISreplaced_genomes (minimal IS extension and merging overlapping IS only) are plotted. Only unitigs with unqiue genome position information (by rounding off to the nearest multiple of 1000) are kept for plotting (as shown in Efaecium32genomes_ext100merge3_1swap_unitig_outdir/kmers_noN/*kmer4plot.txt files). Plots for unitigs can be found in folder Efaecium32genomes_ext100merge3_1swap_unitig_outdir/kmers_noN.

Plot of unitigs that showe rearrangements significantly associated with structure phenotype.

![myNoNintactk_rev0fwd1_unitigs](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/313eba86-1052-48bb-ac29-c4e070cbf11b)

![myNoNintactk_rev1fwd0_unitigs](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/231b818a-df08-46f8-b64b-dba226a23494)

(Above: 191 intact kmers that are in forward orientation in majority of structure "1" genomes, as well as in reverse orientation in majority of structure "0" genomes;
Below: 370 intact kmers that are in reverse orientation in majority of structure "1" genomes, as well as in forward orientation in majority of structure "0" genomes)


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

Genomes asemblies from which genome rearrangements are detected are re-orientated by a chosen gene, i.e. gidA. The location and orientation of gidA in the genomens are obtained by blasting it with multifasta file of genome assemblies.
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
The output file name for the genomes with same orientation is "fixed_genomes.fasta".

Location of IS elements in the genomes are obtained by blasting with target IS sequences. Sequences of more than one target IS element can be placed in the same multifasta file for obtaining their genome locations in all genomes at once.
```
blastn -query example_data/IS_NZ_CP025371.1.fasta \
-subject fixed_genomes.fasta \
-outfmt 6 -out PRN_468_blastIS_out.txt
```

Here, sequences extending 7000bp to both directions from each IS are replaced. IS elements that are no more than 200bp apart (after extension) in each genome are also "merged". Then, each of these "extended and merged" IS element are replaced with shorter placedholder sequences (N x 15). A seperate set of IS-replaced genomes are also produced by enabling performing minimal IS extension (i.e. 100bp) and merging overlapping IS only (i.e. IS that are less than 3 bp apart) through passing string argument "on" to the -s flag.
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

Generate number of unique patterns and p value significance threshold information (For both ext7000_merge200_ISreplaced_genomes set and ext100_merge3_ISreplaced_genomes set):
```
#Run inside corresponding *_ISreplaced_genomes directory

python3 ../scripts/count_patterns.py kmer_patterns_covariate.txt > count_pattern.txt

#count_patterns.py is a script from pyseer package for calculating p-value threshold using Bonferroni correction. 

```
Extract kmers with p value below the significance threshold (For both ext7000_merge200_ISreplaced_genomes set and ext100_merge3_ISreplaced_genomes set):
```
#Run inside corresponding *_ISreplaced_genomes directory

#For ext7000_merge200_ISreplaced_genomes set
awk '{ if ($4 <= 1.53E-04) { print } }' PRN_468_ext7000merge200_k200_MAF0.05_nopopctrl > sigk_pyseer.txt

#For ext100_merge3_ISreplaced_genomes set
awk '{ if ($4 <=  1.62E-04) { print } }' PRN_468_ext100merge3_k200_MAF0.05_nopopctrl > sigk_pyseer.txt
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

Then, these kmers are blasted with the original genome set for studying potential genome rearrangements that are captured by them, implemented by the following script:

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

2) Plotting intact kmers without N for visualising sequence content of rearrangement :

Genome position of intact kmers without N from /ext100_merge3_ISreplaced_genomes (merging overlapping IS only) are plotted. Only kmers with unqiue genome position information (by rounding off to the nearest multiple of 1000) are used in the plot. Plots for intact kmers can be found in folder PRN_468_ext100_merge3_outdir/kmers_noN.

Plots of intact kmers that show sequence rearrangements that are significantly associated with PRN expression phenotype.

<img width="1053" alt="Screenshot 2023-12-13 at 01 34 22" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/cbcc34a4-fc21-442a-b9ca-6849f483db66">

(Left: 7 intact kmers that are in forward orientation in half of PRN+ (1) genomes, as well as reverse in majority of PRN- (0) genomes;
Right: 5 intact kmers that are in reverse orientation in half of PRN+ (1) genomes, as well as forward in majority of PRN- (0) genomes)



**Important notes:**

Some of the significant intack kmers without N contain sequence of pertactin autotransporter (indicated by black arrows). These kmers are not found using ext7000_merge200_ISreplaced_genomes. This is because the gene pertactin autotransporter is located immediately next to an IS element in _pertusis_ genomes, and any genome rearrangement that sits completely within the "replaced IS" region will not be detected. 

Ref: Lefrancq, N., Bouchez, V., Fernandes, N., Barkoff, A.M., Bosch, T., Dalby, T., Åkerlund, T., Darenberg, J., Fabianova, K., Vestrheim, D.F. and Fry, N.K., 2022. Global spatial dynamics and vaccine-induced fitness changes of Bordetella pertussis. Science Translational Medicine, 14(642), p.eabn3253.

