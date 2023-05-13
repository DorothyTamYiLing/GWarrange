Click [here](https://drive.google.com/drive/folders/1RQU1a7kxcVSsIQr94MwSKgRf9PC0YnZO?usp=share_link) for downloading genomes files required for the tutorials. After installing genome_rearangement with git clone, place the downloaded genomes files in ~/genome_rearrangement/example_data for later use.

# Tutotrial 1

This tutorial is based a subset of _Bordetella pertussis_ genomes as described in Weigand _et al_. 2019), in which chromosome structures were defined by exhaustive pairwise alignment. A subset of 47 genomes that display two different chromosome structures (18 genomes with structure 1 and 29 genomes with structure 2) (See Fig. 1) were used. Structure phenotype of two pairs of genomes were swapped for demonstration purpose.  

<img width="1159" alt="Screenshot 2023-04-02 at 10 32 23 PM" src="https://user-images.githubusercontent.com/34043893/229380077-f0c15dba-7ed3-4fc3-a0d7-6d9a5f52b556.png">
Fig 1: Two different chromosome structures were found among 47 _Bordetella pertussis_ genomes. 

First, genomes asemblies from which detecting genome rearrangements are detected are re-orientated by a chosen gene, in the case of Boredetella pertussis is gidA since it is the first gene after origin of replication. The location and orientation of gidA in the genomens are obtained by blasting it with multifasta file of genome assemblies.
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

Genome rearrangments in _Bordetella pertussis_ are belived to be largely mediated by homologous recombination of insertion sequence (IS) elements (such as IS481 and IS110). Location of IS elements in the genomes are obtained by blasting. Sequences of more than one IS elements can be put in the same multifasta file for obtaining genome locations for all at once.
```
blastn -query ./example_data/IS_NZ_CP025371.1.fasta \
-subject fixed_genomes.fasta \
-outfmt 6 -out blastIS_out.txt
```

In addition, genome rearrangments in _Bordetella pertussis_ have also been observed to be mediated by homologous recombination of sequence blocks that consist of one or more IS elements. These duplicated sequence blocks are found throughout the genome and can be as large as several thousand base pairs in size. To ensure sensitivity in detecting genome rearrangements, it is advised to replace the these sequence blocks **completely** with placeholder sequence. Without additional informaiton of the actual size of the homologous sequence blocks, sequences extending several thousands base pairs to both direction from each IS can be replaced.

Here, sequences extending 7000bp to both direction from each IS were replaced. IS elements that were no more than 200bp apart (after extension) in each genome were also "merged". Then, each of these "extended and merged" IS element were replaced with shorter placedholder sequences (N x 15). A seperate set of IS-replaced genomes were also produced by enabling performing minimal extension (i.e. 100bp) and IS replacement with merging overlapping IS only (i.e. IS that are less than 3 bp apart) through passing string argument "on" to the -s flag.
```
bash ./scripts/merge_replace_IS.sh -g fixed_genomes.fasta -i blastIS_out.txt -e 7000 -m 200 -s "on"
```

Each set of IS-replaced genomes using different IS merging and extending parameters were output into new directories "ext7000_merge200_ISreplaced_genomes" and "ext100_merge3_ISreplaced_genomes" (merging overlapping IS only) respectively.

Statistics on the sizes and distance between each pair of adjacent "merged" IS are printed in files *_mergedISstat.txt.

Prior to GWAS, each set of IS-replaced genomes using different IS merging and extending parameters were used for generating kmers. Here, we use the kmer generation tool fsm-lite to generate kmers from genomes in directory /ext7000_merge200_ISreplaced_genomes.

```
cd ./ext7000_merge200_ISreplaced_genomes

#generating fsm-ite input file
for f in *_ext7000_merge200_ISreplaced.fasta; do id=$(basename "$f" _ext7000_merge200_ISreplaced.fasta); echo $id $f; done > clus1clus2_47_input.list

#generating kmers with size of 200 bases with minor allel frequency 0.05
fsm-lite -l clus1clus2_47_input.list -v -s 3 -S 44 -t tmp -m 200 -M 200 | gzip - > clus1clus2_47_ext7000merge200_k200_output.txt.gz
```

Then, a kmer-based GWAS was conducted using pyseer with an aim to identify kmers whose presence-absence patterns are associated with chromosome structures (_i.e._ phenotype). 

```
#adding header to phenotype file for pyseer input format
echo "samples binary" | cat - ../example_data/clus1clus2_pheno.txt > ../example_data/clus1clus2_pheno_4pyseer.txt

#run pyseer
pyseer --phenotypes ../example_data/clus1clus2_pheno_4pyseer.txt \
--kmers clus1clus2_47_ext7000merge200_k200_output.txt.gz \
--no-distances \
--min-af 0.05 --max-af 0.95 \
--print-samples --output-patterns kmer_patterns.txt \
> clus1clus2_47_ext7000merge200_k200_MAF0.05_nopopctrl
```

Generate number of unique pattterns and p value significance threshold information:
```
/path/to/scripts/count_patterns.py kmer_patterns.txt > count_pattern.txt
```
Extract kmers with p value below the the significance threshod:
```
awk '{ if ($4 <= 5.21E-04) { print } }' clus1clus2_47_ext7000merge200_k200_MAF0.05_nopopctrl > sigk_pyseer.txt
```

26,665 kmers were found to be significantly associated with the structural phenotype. The sequences of the kmers were extracted and placed in a multifasta file.

Extract significant kmer sequences and convert then into fasta format
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
#run in the fist level of /genome_rearrangement directory

bash ./scripts/main.sh -k ./ext100_merge7000_ISreplaced_genomes/sigk_seq.fasta \
-g ./example_data/clus1clus2_47.fna \
-p ./example_data/clus1clus2_pheno.txt -l 200 -d 10000 -f 30 \
-o ./clus1clus2_47_ext100_merge7000_outdir -s 4300 -x 2 -y 1000

```

Note that the value used for -d parameter should be larger than the "Maximum size of merged ISs" value in "ext7000_merge200_mergedISstat.txt".

**Visualising genome rearrangements that are captured by kmer**

756 kmers were found to be split (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes (in mysplitk_out.txt).

1) Plotting split kmers for visualising rearrangement boundaries

Since kmers contain highly redundant information, only kmers with unique information (genome position, case and control count and proportion) were kept. They can be found in output file myshort_splitk_out_uniq.txt.

Four rearrangement boundaries were found, and they potentially refer to two inversion events, i.e. between 43000bp and 3600000bp, as well as between 1500000bp and 2500000bp, one inversion nested within the other. The four boundaries can be indicated by eight different significant split kmers that were mapped to each of the boundary and split in case/control genomes (plots of four split kmers were shown below as examples). Full information of these kmers can be found in output file mysplitk_out.txt.

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

Genome position of intact kmers without N from /ext100_merge3_ISreplaced_genomes (merging overlapping IS only) were plotted. Only kmers with unqiue genome position information (by rounding off to the nearest multiple of 1000) 

Plot of intact kmers that show rearrangements in two genome regions that are significantly associated with the structural phenotype.

![myNoNintactk_rev0fwd1](https://user-images.githubusercontent.com/34043893/231439466-2dce018c-89f4-4873-a4a2-486f655af203.png)

![myNoNintactk_rev1fwd0](https://user-images.githubusercontent.com/34043893/231439513-45c179bf-8bc6-4d1b-98b6-fe8f896ce57e.png)

Height of arrows corresponds to proporiton of case/control genomes.

(Above: 15 intact kmers that are in forward orientation in majority of structure "1" genomes, as well as reverse in majority of structure "0" genomes;
Below: 17 intact kmers that are in reverse orientation in majority of structure "1" genomes, as well as forward in majority of structure "0" genomes)

Colour indices refer to the "colour index" column in the corresponding *kmer4plot.txt file (with the same prefix), hence the corresponding kmers.

Reference: Weigand, M.R., Williams, M.M., Peng, Y., Kania, D., Pawloski, L.C., Tondella, M.L. and CDC Pertussis Working Group, 2019. Genomic survey of Bordetella pertussis diversity, United States, 2000–2013. Emerging infectious diseases, 25(4), p.780.


# Tutorial 2

This tutorial is based on 468 _Bordetella pertussis_ genomes with pertactin (PRN) expression information. Among them, 165 genomes show the presence of pertactin expression and 303 show absence. 

First, genomes asemblies from which detecting genome rearrangements are detected are re-orientated by a chosen gene, i.e. gidA. The location and orientation of gidA in the genomens are obtained by blasting it with multifasta file of genome assemblies.
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

Here, sequences extending 7000bp to both direction from each IS were replaced. IS elements that were no more than 200bp apart (after extension) in each genome were also "merged". Then, each of these "extended and merged" IS element were replaced with shorter placedholder sequences (N x 15). A seperate set of IS-replaced genomes were also produced by enabling performing minimal IS extension (i.e. 100bp) and IS replacement with merging overlapping IS only (i.e. IS that are less than 3 bp apart) through passing string argument "on" to the -s flag.
```
bash ./scripts/merge_replace_IS.sh -g fixed_genomes.fasta -i blastIS_out.txt -e 7000 -m 200 -s "on"
```
Prior to GWAS, each set of IS-replaced genomes using different IS merging and extending parameters were used for generating kmers. Here, we use the kmer generation tool fsm-lite to generate kmers from genomes in directory /ext7000_merge200_ISreplaced_genomes.
```
cd ./ext7000_merge200_ISreplaced_genomes
```

#generating fsm-ite input file
```
for f in *_ext100_merge3_ISreplaced.fasta; do id=$(basename "$f" _ext100_merge3_ISreplaced.fasta); echo $id $f; done > PRN_468_input.list
```

#generating kmers with size of 200 bases with minor allel frequency 0.05
```
fsm-lite -l PRN_468_input.list -v -s 24 -S 444 -t tmp -m 200 -M 200 | gzip - > PRN_468_ext100merge3_k200_output.txt.gz
```
Then, a kmer-based GWAS was conducted using pyseer with an aim to identify kmers whose presence-absence patterns are associated with PRN expression phenotype. Population structure is controlled by phylogenetic similarity matrix.

```
#adding header to phenotype file for pyseer input format
echo "samples binary" | cat - ../example_data/prn_status_pheno_4pyseer.txt > ../example_data/prn_status_pheno_4pyseer.txt

#run pyseer
pyseer --lmm --phenotypes ../example_data/prn_status_pheno_4pyseer.txt \
--kmers PRN_468_ext7000merge200_k200_output.txt.gz \
--similarity ../example_data/PRN_USA_468_ClfML_tree_similarity.tsv \
--min-af 0.05 --max-af 0.95 \
--print-samples --output-patterns kmer_patterns.txt \
> PRN_468_ext7000merge200_k200_MAF0.05_nopopctrl
```

Generate number of unique pattterns and p value significance threshold information:
```
/path/to/scripts/count_patterns.py kmer_patterns.txt > count_pattern.txt
```
Extract kmers with p value below the the significance threshod:
```
awk '{ if ($4 <= 1.53E-04) { print } }' PRN_468_ext7000merge200_k200_MAF0.05_nopopctrl > sigk_pyseer.txt
```

18,626 kmers were found to be significantly associated with the PRN expression. The sequences of the kmers were extracted and placed in a multifasta file.

Extract significant kmer sequences and convert then into fasta format
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
#run in the fist level of /genome_rearrangement directory

bash ./scripts/main.sh -k ./ext7000_merge200_ISreplaced_genomes/sigk_seq.fasta \
-g ./example_data/PRN_468.fna \
-p ./example_data/prn_status_pheno.txt -l 200 -d 5000 -f 30 \
-o ./PRN_468_ext7000_merge200_outdir -s 110000 -x 2 -y 1000

```

Note that the value used for -d parameter should be larger than the "Maximum size of merged ISs" value in "ext7000_merge200_mergedISstat.txt".

**Visualising genome rearrangements that are captured by kmer**

126 kmers were found to be intact in majority of case genomes and split in half of control genomes when mapped to the original genomes. They refer to one rearrangement boundary that potentially belongs to an inversion event between 1800000bp and 2200000bp in the genomes. Full information of these kmers can be found in output file mysplitk_out.txt. One of such split kmers is shown below:

![kmer970_plot-1](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/f9190541-f0df-4fbf-b645-ede2c26bdbab)

Height of arrows corresponds to proportion of case/control genomes.

Genome position of intact kmers without N from /ext100_merge3_ISreplaced_genomes (merging overlapping IS only) were plotted. Only kmers with unqiue genome position information (by rounding off to the nearest multiple of 1000) 

Plot of intact kmers that show rearrangements in two genome regions that are significantly associated with the structural phenotype.

![myNoNintactk_rev0fwd1](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/d0022eff-2cec-4941-9ba4-c8fc69567924)

![myNoNintactk_rev1fwd0](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/b39e83fa-a961-4d7a-8de7-fa8f04e4b32d)

Height of arrows corresponds to proporiton of case/control genomes.

(Above: 4 intact kmers that are in forward orientation in majority of structure "1" genomes, as well as reverse in majority of structure "0" genomes;
Below: 3 intact kmers that are in reverse orientation in majority of structure "1" genomes, as well as forward in majority of structure "0" genomes)

Colour indices refer to the "colour index" column in the corresponding *kmer4plot.txt file (with the same prefix), hence the corresponding kmers.

