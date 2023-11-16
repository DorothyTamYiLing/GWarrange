
<img width="591" alt="Screenshot 2023-09-18 215805" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/f3501563-7c65-4c04-98c5-12a13ed3c3ad">

Genome rearrangement pipeline summary chart

Click [here](https://drive.google.com/drive/folders/1RQU1a7kxcVSsIQr94MwSKgRf9PC0YnZO?usp=share_link) for downloading genomes files required for the tutorials. After installing genome_rearangement with git clone, place the downloaded genomes files in ~/genome_rearrangement/example_data for later use.

# Tutorial 1

This tutorial is based a subset of _Bordetella pertussis_ genomes as described in Weigand _et al_. 2019, in which chromosome structures were defined by MAUVE exhaustive pairwise alignment. A subset of 47 genomes that display two different chromosome structures (18 genomes with structure 1 and 29 genomes with structure 2) (See Fig. 1) were used in GWAS with structure information as phenotype. Structure phenotype of two pairs of genomes were swapped for demonstration purpose.  

<img width="1159" alt="Screenshot 2023-04-02 at 10 32 23 PM" src="https://user-images.githubusercontent.com/34043893/229380077-f0c15dba-7ed3-4fc3-a0d7-6d9a5f52b556.png">
Fig 1: Two different chromosome structures were found among 47 _Bordetella pertussis_ genomes. 

##

First, the genome assemblies multifasta file is prepared by concatenating genome fasta files.
```
#concatenating genome fasta files for use
cd ./example_data/example_genomes/clus1clus2_47_genomes
cat *fasta.gz > ../../clus1clus2_47.fna.gz
```

Genomes asemblies from which genome rearrangements are detected are re-orientated by a chosen gene. In the case of Boredetella pertussis, the gene is gidA since it is the first gene after origin of replication. The location and orientation of gidA in the genomes are obtained by blasting it with multifasta file of genome assemblies.
```
#unzip the genome file if neccesasry
gunzip ./example_data/clus1clus2_47.fna.gz

#blast gidA with genomes
blastn -query ./example_data/gidA.fasta \
-subject ./example_data/clus1clus2_47.fna \
-outfmt 6 -out clus1clus2_47_gidA_out.txt
```

Then, genome assemblies are re-orientated according to the position and orientation of gidA in the genomes, using the script fix_genome.py:
```
python3 ./scripts/fix_genome.py --input ./example_data/clus1clus2_47.fna --mycoor clus1clus2_47_gidA_out.txt
```

The output file name for the genomes with same orientation is "fixed_genomes.fasta".

Genome rearrangments in _Bordetella pertussis_ are belived to be largely mediated by homologous recombination of insertion sequence (IS) elements (such as IS481 and IS110). Location of IS elements in the genomes are obtained by blasting with target IS sequences. Sequences of more than one target IS elements can be placed in the same multifasta file for obtaining their genome locations in all genomes at once.
```
blastn -query ./example_data/IS_NZ_CP025371.1.fasta \
-subject fixed_genomes.fasta \
-outfmt 6 -out clus1clus2_47_blastIS_out.txt
```

In addition, genome rearrangments in _Bordetella pertussis_ have also been observed to be mediated by homologous recombination of sequence blocks that consist of one or more IS elements. These duplicated sequence blocks are found throughout the genome and can be as large as several thousand base pairs in size. To ensure sensitivity in genome rearrangements detection, it is advised to replace these sequence blocks **completely** with shorter placeholder sequence. Without additional information of the actual size of the homologous sequence blocks, sequences extending several thousands base pairs to both direction from each IS can be replaced.

Here, sequences extending 7000bp to both direction from each IS were replaced. IS elements that were no more than 200bp apart (after extension) in each genome were also "merged". Then, each of these "extended and merged" IS elements were replaced with shorter placedholder sequences (N x 15). A seperate set of IS-replaced genomes were also produced by enabling performing minimal IS extension (i.e. 100bp to both direction) and merging overlapping IS only (i.e. ISs that are less than 3bp apart) through passing string argument "on" to the -s flag.
```
bash ./scripts/merge_replace_IS.sh -g fixed_genomes.fasta -i clus1clus2_47_blastIS_out.txt -e 7000 -m 200 -s "on"
```

Each set of IS-replaced genomes using different IS merging and extending parameters were output into new directories "ext7000_merge200_ISreplaced_genomes" and "ext100_merge3_ISreplaced_genomes" respectively.

Statistics on the sizes and distances between each pair of adjacent "extended and merged" IS are printed in files *_mergedISstat.txt.

Prior to GWAS, each set of IS-replaced genomes using different IS merging and extending parameters were used for generating kmers.

```
#For ext7000_merge200_ISreplaced_genomes set
cd ./ext7000_merge200_ISreplaced_genomes

#generating fsm-ite input file
for f in *_ext7000_merge200_ISreplaced.fasta; do id=$(basename "$f" _ext7000_merge200_ISreplaced.fasta); echo $id $f; done > clus1clus2_47_input.list

#generating kmers with size of 200 bases with minor allel frequency 0.05
fsm-lite -l clus1clus2_47_input.list -v -s 3 -S 44 -t tmp -m 200 -M 200 | gzip - > clus1clus2_47_ext7000merge200_k200_output.txt.gz

########################################################################

#For ext100_merge3_ISreplaced_genomes set
cd ./ext100_merge3_ISreplaced_genomes

for f in *_ext100_merge3_ISreplaced.fasta; do id=$(basename "$f" _ext100_merge3_ISreplaced.fasta); echo $id $f; done > clus1clus2_47_input.list

#generating kmers with size of 200 bases with minor allel frequency 0.05
fsm-lite -l clus1clus2_47_input.list -v -s 3 -S 44 -t tmp -m 200 -M 200 | gzip - > clus1clus2_47_ext100merge3_k200_output.txt.gz

```

Then, a kmer-based GWAS was conducted using pyseer with an aim to identify kmers whose presence-absence patterns were associated with chromosome structures phenotype. Population structure is not contolled.

```
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

Generate number of unique pattterns and p value significance threshold information:
```
#count_patterns.py is a script from pyseer package for calculating p-value threshold using Bonferroni correction. To access pyseer scripts, one needs to have cloned the pyseer github repository and go to scripts/directory.
../scripts/count_patterns.py kmer_patterns.txt > count_pattern.txt
```
Extract kmers with p value below the the significance threshold:
```
#For ext7000_merge200_ISreplaced_genomes set
awk '{ if ($4 <= 4.59E-04) { print } }' clus1clus2_47_ext7000merge200_k200_MAF0.05_nopopctrl > sigk_pyseer.txt

#For ext100_merge3_ISreplaced_genomes set
awk '{ if ($4 <= 5.62E-04) { print } }' clus1clus2_47_ext100merge3_k200_MAF0.05_nopopctrl > sigk_pyseer.txt
```

The sequences of kmers that were found to be significantly associated with the structural phenotyp were extracted and placed in a multifasta file.

Extract significant kmer sequences and convert them into fasta format:
```
#get the seqeunce only
awk '{print $1}' sigk_pyseer.txt > sigk_seq.txt 

#create multifasta file for significant kmer sequences
number=$(cat sigk_seq.txt | wc -l)

rm header.txt   #remove any existing header file

START=1
let "END=$number" 
 
for (( c=$START; c<=$END; c++ ))
do
	echo ">kmer""$c " >> header.txt
done

paste -d \\n header.txt sigk_seq.txt > sigk_seq.fasta
```

Then, these kmers were blasted with the original genome set for studying potential genome rearrangment that are captured by them, implemented by the following script:

```
#run in the first level of /genome_rearrangement directory

#For ext7000_merge200_ISreplaced_genomes set
bash ./scripts/main.sh -k ./ext7000_merge200_ISreplaced_genomes/sigk_seq.fasta \
-g ./example_data/clus1clus2_47.fna.gz \
-p ./example_data/clus1clus2_pheno.txt -d 110000 -f 30 \
-o ./clus1clus2_47_ext7000_merge200_outdir -s 4300 -x 2 -y 1000

#For ext100_merge3_ISreplaced_genomes set
bash ./scripts/main.sh -k ./ext100_merge3_ISreplaced_genomes/sigk_seq.fasta \
-g ./example_data/clus1clus2_47.fna \
-p ./example_data/clus1clus2_pheno.txt -d 5000 -f 30 \
-o ./clus1clus2_47_ext100_merge3_outdir -s 4300 -x 2 -y 1000

```

Note that the value used for -d parameter should be larger than the "Maximum size of merged ISs" value in the corresponding *mergedISstat.txt file.

**Visualising genome rearrangements that are captured by kmer**

From genome set with 7000bp extension and 200bp merging, 978 kmers were found to be split (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes.

1) Plotting split kmers for visualising rearrangement boundaries

Since kmers contain highly redundant information, only kmers with unique information (genome position, case and control count and proportion) were kept. They can be found in output file myshort_splitk_out_uniq.txt.

Four rearrangement boundaries were found, and they potentially refer to two inversion events, i.e. between 43000bp and 3600000bp, as well as between 1500000bp and 2500000bp, one inversion nested within the other. The four boundaries can be indicated by sixteen different significant split kmers that were mapped to each of the boundary, split in case/control genomes, and in forward/reverse orientation (plots of four split kmers were shown below as examples). Full information of these kmers can be found in output file mysplitk_out.txt.

Inversion within genome region 43000 and 3600000, 43000bp boundary, kmer being intact in case genomes and split in control genomes, in forward orientation:

<img width="927" alt="Screenshot 2023-05-07 at 9 37 35 PM" src="https://user-images.githubusercontent.com/34043893/236681154-15c6b72f-9c61-40d9-9fd4-04ae35d22620.png">

Inversion within genome region 43000 and 3600000, 3600000bp boundary, kmer being intact in control genomes and split in case genomes, in reverse orientation:

<img width="910" alt="Screenshot 2023-05-07 at 9 45 08 PM" src="https://user-images.githubusercontent.com/34043893/236681643-a20c9927-79fb-47b3-8a8b-ebd83a239528.png">

Inversion within genome region 1500000 and 2500000, 1500000bp boundary, kmer being intact in control genomes and split in case genomes, forward orientation kmer:

<img width="906" alt="Screenshot 2023-05-07 at 9 50 43 PM" src="https://user-images.githubusercontent.com/34043893/236681673-c8ac9e6a-fa77-4314-af17-022e6b851552.png">

Inversion within genome region 1500000 and 2500000, 2500000bp boundary, kmer being intact in case genomes and split in control genomes, reverse orientation kmer:

<img width="904" alt="Screenshot 2023-05-07 at 9 48 39 PM" src="https://user-images.githubusercontent.com/34043893/236681684-c2a4eaf5-0460-4243-8d3c-f87988834efa.png">

Height of arrows corresponds to proporiton of case/control genomes.

2) Plotting intact kmers without N for visualising sequence content of rearrangement :

Genome position of intact kmers without N from /ext100_merge3_ISreplaced_genomes (minimal IS extension and merging overlapping IS only) were plotted. Only kmers with unqiue genome position information (by rounding off to the nearest multiple of 1000) were kept for plotting (as shown in *kmer4plot.txt files). 

Plot of intact kmers that showed rearrangements in two genome regions that were significantly associated with structural phenotype.

![myNoNintactk_rev0fwd1](https://user-images.githubusercontent.com/34043893/231439466-2dce018c-89f4-4873-a4a2-486f655af203.png)

![myNoNintactk_rev1fwd0](https://user-images.githubusercontent.com/34043893/231439513-45c179bf-8bc6-4d1b-98b6-fe8f896ce57e.png)

Height of arrows corresponds to proporiton of case/control genomes.

(Above: 15 intact kmers that are in forward orientation in majority of structure "1" genomes, as well as in reverse orientation in majority of structure "0" genomes;
Below: 17 intact kmers that are in reverse orientation in majority of structure "1" genomes, as well as in forward orientation in majority of structure "0" genomes)

Colour indices refer to the "colour index" column in the corresponding *kmer4plot.txt file (with the same prefix), hence the corresponding kmers.

**Important notes:**

Inversion boundaries between genome region 1500000 and 2500000 were not detected using ext100_merge3_ISreplaced_genomes. This is because the whole homologous sequence block/IS clusters at these boundaries are not completely replaced by shorter placeholder sequences, this leads to the presence of homologous sequence within flanking sequences, hence reducing senesitivity in detecting rearrangement boundaries. To ensure replacing the whole homologous sequence block/IS clusters by short placeholder sequences, we can extending the genome coordinates of each repeated sequence for a number of base pairs in both directions, and/or to merge repeated sequences that are less than a number of base pairs distance apart.

Reference: Weigand, M.R., Williams, M.M., Peng, Y., Kania, D., Pawloski, L.C., Tondella, M.L. and CDC Pertussis Working Group, 2019. Genomic survey of Bordetella pertussis diversity, United States, 2000–2013. Emerging infectious diseases, 25(4), p.780.


# Tutorial 2

_Enterococcus faecium_'s genomes are known to be enriched with IS elements, which could play important role in their genome structure's diversification (Leavis et al. 2007). Genome structure of 75 Enterococcus faecium were characterised by socru (Page et al. 2020), Among which, a subset of 32 genomes display two different chromosome structures (21 genomes with structure "0" and 11 genomes with structure "1") (See Fig. ) were used in GWAS, with structure information as phenotype. Structure phenotype of two pairs of genomes were swapped for demonstration purpose. 

![Screenshot 2023-08-01 222347](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/0f0c348f-936b-4983-be54-180cc6b8d838)
Fig : Two different chromosome structures were found among 32 _Enterococcus faecium_ genomes. 

First, genomes asemblies from which genome rearrangements are detected are re-orientated by a chosen gene, i.e. dnaA. The location and orientation of dnaA in the genomens are obtained by blasting it with multifasta file of genome assemblies.

```
#unzip the genome file if neccesasry
gunzip ./example_data/32genomes.fna.gz

#blast dnaA with genomes
blastn -query ./example_data/dnaA.fasta \
-subject ./example_data/32genomes.fna \
-outfmt 6 -out 32genomes_dnaA_out.txt
```

Then, genome assemblies are re-orientated according to the position and orientation of dnaA in the genomes, using the script fix_genome.py:

```
python3 ./scripts/fix_genome.py --input ./example_data/32genomes.fna --mycoor 32genomes_dnaA_out.txt
```
The output file name for the genomes with same orientation is "fixed_genomes.fasta".

Location of IS elements in the genomes are obtained by blasting with target IS sequences. Sequences of more than one target IS elements can be placed in the same multifasta file for obtaining their genome locations in all genomes at once.

```
blastn -query ./example_data/IS30_IS1252_in_HOU503_657692to658645.fasta \
-subject fixed_genomes.fasta \
-outfmt 6 -out blastIS30_IS1252_in_HOU503_32genomes_out.txt
```

Here, sequences extending 7000bp to both direction from each IS were replaced. IS elements that were no more than 200bp apart (after extension) in each genome were also "merged". Then, each of these "extended and merged" IS element were replaced with shorter placedholder sequences (N x 15). A seperate set of IS-replaced genomes were also produced by enabling performing minimal IS extension (i.e. 100bp) and merging overlapping IS only (i.e. IS that are less than 3 bp apart) through passing string argument "on" to the -s flag.
```
bash ./scripts/merge_replace_IS.sh -g fixed_genomes.fasta -i blastIS30_IS1252_in_HOU503_32genomes_out.txt -e 7000 -m 200 -s "on"
```
Prior to GWAS, each set of IS-replaced genomes using different IS merging and extending parameters were used for generating kmers and unitigs. Here, we use the kmer generation tool fsm-lite and unitig generation tool unitig-caller to generate kmers and unitigs from genomes in directory /ext100_merge3_ISreplaced_genomes.
```
cd ./ext100_merge3_ISreplaced_genomes

#generating fsm-ite input file
for f in *_ext100_merge3_ISreplaced.fasta; do id=$(basename "$f" _ext100_merge3_ISreplaced.fasta); echo $id $f; done > 32genomes_input.list

#generating kmers with size of 200 bases kmers present in at least 20 samples

fsm-lite -l input.list -v -s 20 -S 30 -t tmp -m 200 -M 200 | gzip - > ext100merge3_k200_output.txt.gz

#generating input file for unitig-caller
ls -d -1 $PWD/*.fasta > input.txt

#running unitig-caller
unitig-caller --call --rtab --pyseer --refs input.txt --out unitigcall_out
```
Then, kmer-based and unitig-based GWAS were conducted using pyseer with an aim to identify kmers/unitigs whose presence-absence patterns were associated with chromosome structures phenotype. Population structure is not contolled.
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
Generate number of unique pattterns and p value significance threshold information:
```
#count_patterns.py is a script from pyseer package for calculating p-value threshold using Bonferroni correction. To access pyseer scripts, one needs to have cloned the pyseer github repository and go to scripts/directory.
../scripts/count_patterns.py kmer_patterns.txt > count_pattern.txt

../scripts/count_patterns.py unitigs_kmer_patterns.txt > count_pattern.txt
```
Extract kmers/unitigs with p value below the the significance threshod:
```
awk '{ if ($4 <= 1.92E-05) { print } }' ext100merge3_k200_min20samp_nopopctrl > sigk_pyseer.txt

awk '{ if ($4 <= 7.62E-06) { print } }' ext100merge3_k200_min20samp_unitigs_nopopctrl > siguni_pyseer.txt
```
448,330 kmers and 3,737 unitigs were found to be significantly associated with chromosome structure. The sequences of which were extracted and placed in respective multifasta file.

Extract significant kmer sequences and convert them into fasta format
```
#get the seqeunce only
awk '{print $1}' sigk_pyseer.txt > sigk_seq.txt 

#create multifasta file for significant kmer sequences
number=$(cat sigk_seq.txt | wc -l)

rm header.txt   #remove any existing header file

START=1
let "END=$number" 
 
for (( c=$START; c<=$END; c++ ))
do
	echo ">kmer""$c " >> header.txt
done

paste -d \\n header.txt sigk_seq.txt > sigk_seq.fasta
```
#Due to the large number of significant kmers, only the kmers that contain "N" and the first 5000 kmers without "N" were used for structural analysis.
```
#classifying kmers into those containing "N" and those do not contain "N"
python3 ../scripts/class_k.py --input sigk_seq.fasta --outdir .

#extract first 5000 kmers without "N"
head -25000 sigk_noN.fasta > sigk_noN_5000.fasta

#combine the kmers with "N" and the first 5000 kmers without "N"
cat sigk_withN.fasta sigk_noN_5000.fasta > sigkwithN_noN5000.fasta
```

Extract significant unitig sequences and convert then into fasta format
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

Then, these kmers and unitigs were blasted with the original genome set for studying potential genome rearrangment that are captured by them, implemented by the following script:

#run in the first level of /genome_rearrangement directory
```
#For ext100_merge3_ISreplaced_genomes set k-mers with N 
bash ./scripts/main.sh -k ./ext100_merge3_ISreplaced_genomes/sigkwithN_noN5000.fasta \
-g ./example_data/Efaecium_32genomes.fna.gz \
-p ./example_data/Efaecium32genomes_pheno_1swap.txt -d 3000 -f 30 \
-o ./Efaecium32genomes_ext100merge3_1swap_withNnoN5000_outdir -s 3000 -x 2 -y 1000

#For ext100_merge3_ISreplaced_genomes set unitigs
bash ./scripts/main.sh -k ./ext100_merge3_ISreplaced_genomes/siguni_seq.fasta \
-g ./example_data/Efaecium_32genomes.fna.gz \
-p ./example_data/Efaecium32genomes_pheno_1swap.txt -d 3000 -f 30 \
-o ./Efaecium32genomes_ext100merge3_1swap_unitig_outdir -s 3000 -x 2 -y 1000
```
Note that the value used for -d parameter should be larger than the "Maximum size of merged ISs" value in "ext100_merge3_mergedISstat.txt".

**Visualising genome rearrangements that are captured by kmer**

219 kmers were found to be split (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes.

1) Plotting split kmers for visualising rearrangement boundaries

Since kmers contain highly redundant information, only kmers with unique information (genome position, case and control count and proportion) were kept. They can be found in output file myshort_splitk_out_uniq.txt.

Two rearrangement boundaries were found, and they potentially refer to a single inversion event between 72000bp and 2100000bp. All of the split kmers were found to be split in case geomes and intact in control genomes. This could be explained by the absence of IS elements in the rearrangement boundaries in case genomes, which has been confirmed by manual sequence check. The two boundaries can be indicated by four different significant split kmers that were mapped to each of the boundary and in forward/reverse orientation (plots of one of these kmers are shown below). Full information of these kmers can be found in output file mysplitk_out.txt.

72000bp boundary, intact kmer in forward orientation:

![kmer98323_plot-1](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/68c8e75f-f445-4bec-b9f2-009c30be3e19)

72000bp boundary, intact kmer in reverse orientation:

![kmer97292_plot-1](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/6ef7c181-0dcc-4f84-a4fe-dc611919d511)

2100000bp boundary, intact kmer in forward orientation:

![kmer99909_plot-1](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/98c528e9-3709-4429-a36f-3c5f2425416f)

2100000bp boundary, intact kmer in reverse orientation:

![kmer98897_plot-1](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/b7d2c49f-fa69-4a32-8b1f-5c7dcaff267d)

Height of arrows corresponds to proporiton of case/control genomes.

No significant split unitigs were identified as unitig-caller was not able to generate unitigs containing placeholder sequences.

2) Plotting intact unitigs (without N) for visualising sequence content of rearrangement :

Genome position of unitigs from /ext100_merge3_ISreplaced_genomes (minimal IS extension and merging overlapping IS only) were plotted. Only unitigs with unqiue genome position information (by rounding off to the nearest multiple of 1000) were kept for plotting (as shown in *kmer4plot.txt files). 

Plot of unitigs that showed rearrangements significantly associated with structural phenotype.

![myNoNintactk_rev0fwd1](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/a9ce4a68-3de4-48d1-b2ab-4675cc2ee8a9)

![myNoNintactk_rev1fwd0](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/2799983d-f31f-4bdc-a200-6e4e58909c00)

Height of arrows corresponds to proporiton of case/control genomes.

(Above: 191 intact kmers that are in forward orientation in majority of structure "1" genomes, as well as in reverse orientation in majority of structure "0" genomes;
Below: 370 intact kmers that are in reverse orientation in majority of structure "1" genomes, as well as in forward orientation in majority of structure "0" genomes)

Colour indices refer to the "colour index" column in the corresponding *kmer4plot.txt file (with the same prefix), hence the corresponding unitigs.

Page, A.J., Ainsworth, E.V. and Langridge, G.C., 2020. socru: typing of genome-level order and orientation around ribosomal operons in bacteria. Microbial Genomics, 6(7).


# Tutorial 3

Leavis, H.L., Willems, R.J.L., van Wamel, W.J.B., Schuren, F.H., Caspers, M.P.M. and Bonten, M.J.M., 2007. Insertion sequence–driven diversification creates a globally dispersed emerging multiresistant subspecies of E. faecium. PLoS pathogens, 3(1), p.e7.

This tutorial is based on 468 _Bordetella pertussis_ genomes with pertactin (PRN) expression information. Among them, 165 genomes show the presence of pertactin expression and 303 show absence. Pertactin expression information is taken from the supplementary material summarised in Lefrancq _et al._ 2022.

First, the genome assemblies multifasta file is prepared by concatenating genome fasta files.
```
#concatenating genome fasta files for use
cd ./example_data/example_genomes/PRN_468
cat *fasta.gz > ../../PRN_468.fna.gz
```

Genomes asemblies from which genome rearrangements are detected are re-orientated by a chosen gene, i.e. gidA. The location and orientation of gidA in the genomens are obtained by blasting it with multifasta file of genome assemblies.
```
#unzip the genome file if neccesasry
gunzip ./example_data/PRN_468.fna.gz

#blast gidA with genomes
blastn -query ./example_data/gidA.fasta \
-subject ./example_data/PRN_468.fna \
-outfmt 6 -out PRN_468_gidA_out.txt
```
Then, genome assemblies are re-orientated according to the position and orientation of gidA in the genomes, using the script fix_genome.py:
```
python3 ./scripts/fix_genome.py --input ./example_data/PRN_468.fna --mycoor PRN_468_gidA_out.txt
```
The output file name for the genomes with same orientation is "fixed_genomes.fasta".

Location of IS elements in the genomes are obtained by blasting with target IS sequences. Sequences of more than one target IS elements can be placed in the same multifasta file for obtaining their genome locations in all genomes at once.
```
blastn -query ./example_data/IS_NZ_CP025371.1.fasta \
-subject fixed_genomes.fasta \
-outfmt 6 -out PRN_468_blastIS_out.txt
```

Here, sequences extending 7000bp to both direction from each IS were replaced. IS elements that were no more than 200bp apart (after extension) in each genome were also "merged". Then, each of these "extended and merged" IS element were replaced with shorter placedholder sequences (N x 15). A seperate set of IS-replaced genomes were also produced by enabling performing minimal IS extension (i.e. 100bp) and merging overlapping IS only (i.e. IS that are less than 3 bp apart) through passing string argument "on" to the -s flag.
```
bash ./scripts/merge_replace_IS.sh -g fixed_genomes.fasta -i PRN_468_blastIS_out.txt -e 7000 -m 200 -s "on"
```
Prior to GWAS, each set of IS-replaced genomes using different IS merging and extending parameters were used for generating kmers.
```
#For ext7000_merge200_ISreplaced_genomes set
cd ./ext7000_merge200_ISreplaced_genomes

#generating fsm-ite input file
for f in *_ext7000_merge200_ISreplaced.fasta; do id=$(basename "$f" _ext7000_merge200_ISreplaced.fasta); echo $id $f; done > PRN_468_input.list

#generating kmers with size of 200 bases with minor allel frequency 0.05

fsm-lite -l PRN_468_input.list -v -s 24 -S 444 -t tmp -m 200 -M 200 | gzip - > PRN_468_ext7000merge200_k200_output.txt.gz

##############################################################################

#For ext100_merge3_ISreplaced_genomes set
cd ./ext100_merge3_ISreplaced_genomes

#generating fsm-ite input file
for f in *_ext100_merge3_ISreplaced.fasta; do id=$(basename "$f" _ext100_merge3_ISreplaced.fasta); echo $id $f; done > PRN_468_input.list

#generating kmers with size of 200 bases with minor allel frequency 0.05
fsm-lite -l PRN_468_input.list -v -s 24 -S 444 -t tmp -m 200 -M 200 | gzip - > PRN_468_ext100merge3_k200_output.txt.gz

```
Then, a kmer-based GWAS was conducted using pyseer with an aim to identify kmers whose presence-absence patterns are associated with PRN expression phenotype. Population structure is controlled by phylogenetic similarity matrix.

```
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

Generate number of unique pattterns and p value significance threshold information:
```
#count_patterns.py is a script from pyseer package for calculating p-value threshold using Bonferroni correction. To access pyseer scripts, one needs to have cloned the pyseer github repository and go to scripts/directory.
./scripts/count_patterns.py kmer_patterns.txt > count_pattern.txt
```
Extract kmers with p value below the the significance threshold:
```
#For ext7000_merge200_ISreplaced_genomes set
awk '{ if ($4 <= ) { print } }' PRN_468_ext7000merge200_k200_MAF0.05_nopopctrl > sigk_pyseer.txt

#For ext100_merge3_ISreplaced_genomes set
awk '{ if ($4 <= ) { print } }' PRN_468_ext100merge3_k200_MAF0.05_nopopctrl > sigk_pyseer.txt
```

The sequences of kmers that were found to be significantly associated with the structural phenotyp were extracted and placed in a multifasta file.

Extract significant kmer sequences and convert them into multifasta format:
```
#get the seqeunce only
awk '{print $1}' sigk_pyseer.txt > sigk_seq.txt 

#create multifasta file for significant kmer sequences
number=$(cat sigk_seq.txt | wc -l)

rm header.txt   #remove any existing header file

START=1
let "END=$number" 
 
for (( c=$START; c<=$END; c++ ))
do
	echo ">kmer""$c " >> header.txt
done

paste -d \\n header.txt sigk_seq.txt > sigk_seq.fasta
```

Then, these kmers were blasted with the original genome set for studying potential genome rearrangment that are captured by them, implemented by the following script:

```
#run in the first level of /genome_rearrangement directory

#For ext7000_merge200_ISreplaced_genomes set
bash ./scripts/main.sh -k ./ext7000_merge200_ISreplaced_genomes/sigk_seq.fasta \
-g ./example_data/PRN_468.fna \
-p ./example_data/prn_status_pheno.txt -d 110000 -f 30 \
-o ./PRN_468_ext7000_merge200_outdir -s 4300 -x 2 -y 1000

#For ext100_merge3_ISreplaced_genomes set
bash ./scripts/main.sh -k ./ext100_merge3_ISreplaced_genomes/sigk_seq.fasta \
-g ./example_data/PRN_468.fna.gz \
-p ./example_data/prn_status_pheno.txt -d 5000 -f 30 \
-o ./PRN_468_ext100_merge3_outdir -s 4300 -x 2 -y 1000

```

Note that the value used for -d parameter should be larger than the "Maximum size of merged ISs" value in the corresponding *_mergedISstat.txt file.

**Visualising genome rearrangements that are captured by kmer**

1) Plotting split kmers for visualising rearrangement boundaries

No split kmer that indicated phenotype-associated rearrangement boundary was detected.

2) Plotting intact kmers without N for visualising sequence content of rearrangement :

Genome position of intact kmers without N from /ext100_merge3_ISreplaced_genomes (merging overlapping IS only) were plotted. Only kmers with unqiue genome position information (by rounding off to the nearest multiple of 1000) 

Plot of intact kmers that show sequence rearrangements that are significantly associated with the PRN expression phenotype.

![myNoNintactk_rev0fwd1](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/d0022eff-2cec-4941-9ba4-c8fc69567924)

![myNoNintactk_rev1fwd0](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/b39e83fa-a961-4d7a-8de7-fa8f04e4b32d)

Height of arrows corresponds to proporiton of case/control genomes.

(Above: 4 intact kmers that are in forward orientation in half of PRN+ (1) genomes, as well as reverse in majority of PRN- (0) genomes;
Below: 3 intact kmers that are in reverse orientation in half of PRN+ (1) genomes, as well as forward in majority of PRN- (0) genomes)

Colour indices refer to the "colour index" column in the corresponding *kmer4plot.txt file (with the same prefix), hence the corresponding kmers.

**Important notes:**

Some of the significant intack kmers without N (e.g. kmer2189, kmer1450, kmer2126) contains sequence of pertactin autotransporter. These kmers was not found using ext7000_merge200_ISreplaced_genomes. This is because the gene pertactin autotransporter is located immediately next to an IS element in pertusis genomes, and any genome rearrangement that sits completely within the "replaced IS" region will not be detected. 

Ref: Lefrancq, N., Bouchez, V., Fernandes, N., Barkoff, A.M., Bosch, T., Dalby, T., Åkerlund, T., Darenberg, J., Fabianova, K., Vestrheim, D.F. and Fry, N.K., 2022. Global spatial dynamics and vaccine-induced fitness changes of Bordetella pertussis. Science Translational Medicine, 14(642), p.eabn3253.

