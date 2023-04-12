Click [here] (https://drive.google.com/drive/folders/1RQU1a7kxcVSsIQr94MwSKgRf9PC0YnZO?usp=share_link) for downloading genomes files required for the tutorials. After installing genome_rearangement with git clone, place the downloaded genomes files in ~/genome_rearrangement/example_data for later use.

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

In addition, genome rearrangments in _Bordetella pertussis_ have also been observed to be mediated by homologous recombination of sequence blocks that consist of multiple IS elements, which are usually found at the start and end of the sequence block. These sequence blocks can be as large as several thousand base pair in size. In order to increase the chance for detecting of genome rearrangements, it is advised to replace the these sequence blocks comlpetely with placeholder sequence. 

Here, IS elements that were no more than 7000bp apart in each genome were "merged". Then, each of these "merged" IS element were replaced with shorter placedholder sequences (N x 15). Each pair of IS cordinates (start and end genome position) were extended by 100bp (default) on each side to ensure complete mask of the IS. A separaet set of IS-replaced genomes were also produced by enabling performing IS replacement with merging overlapping IS only (i.e. IS that are less than 3 bp apart) through passing string argument "on" to the -s flag.
```
bash ./scripts/merge_replace_IS.sh -g fixed_genomes.fasta -i blastIS_out.txt -m 7000 -s "on"
```

Each set of IS-replaced genomes using different IS merging and extending parameters were output into new directories "ext100_merge7000_ISreplaced_genomes" and "ext100_merge3_ISreplaced_genomes" (merging overlapping IS only) respectively.

Prior to GWAS, each set of IS-replaced genomes using different IS merging and extending parameters were used for generating kmers. Here, we use the kmer generation tool fsm-lite to generate kmers from genomes in directory /ext100_merge7000_ISreplaced_genomes.

```
cd ./ext100_merge7000_ISreplaced_genomes

#generating fsm-ite input file
for f in *_ext100_merge7000_ISreplaced.fasta; do id=$(basename "$f" _ext100_merge7000_ISreplaced.fasta); echo $id $f; done > clus1clus2_47_input.list

#generating kmers with size of 200 bases with minor allel frequency 0.05
fsm-lite -l clus1clus2_47_input.list -v -s 3 -S 44 -t tmp -m 200 -M 200 | gzip - > clus1clus2_47_ext100merge7000_k200_output.txt.gz
```

Then, a kmer-based GWAS was conducted using pyseer with an aim to identify kmers whose presence-absence patterns are associated with chromosome structures (_i.e._ phenotype). 

```
#adding header to phenotype file for pyseer input format
echo "samples binary" | cat - ../example_data/clus1clus2_pheno.txt > ../example_data/clus1clus2_pheno_4pyseer.txt

#run pyseer
pyseer --phenotypes ../example_data/clus1clus2_pheno_4pyseer.txt \
--kmers clus1clus2_47_ext100merge7000_k200_output.txt.gz \
--no-distances \
--min-af 0.05 --max-af 0.95 \
--print-samples --output-patterns kmer_patterns.txt \
> clus1clus2_47_ext100merge7000_k200_MAF0.05_nopopctrl
```

Generate number of unique pattterns and p value significance threshold information:
```
/path/to/scripts/count_patterns.py kmer_patterns.txt > count_pattern.txt
```
Extract kmers with p value below the the significance threshod:
```
awk '{ if ($4 <= 5.21E-04) { print } }' clus1clus2_47_ext5merge7000_k200_MAF0.05_nopopctrl > sigk_pyseer.txt
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
-p ./example_data/clus1clus2_pheno.txt -l 200 -d 7000 -f 30 \
-o ./clus1clus2_47_ext100_merge7000_outdir -s 4300 -x 2 -y 1000

```

756 kmers were found to be split (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes (in mysplitk_out.txt).


**Visualising genome rearrangements that are captured by kmer**

1) Plotting split kmers for visualising rearrangement boundaries

Since kmers contain highly redundant information, only kmers with unique information (genome position, case and control count and proportion) were kept. They can be found in output file myshort_splitk_out_uniq.txt.

Rearrangement boundaries that potentially refer to two inversion events, one nested within the other, were detected. Here, in total, seven significant split kmers contained full information for the boundaries of two inversion events (plot of two split kmers were shown below as examples). Full information of these kmers can be found in output file mysplitk_out.txt.

Inversion within genome region 43000 and 3600000, 3600000bp boundary, kmer being intact in control genomes and split in case genomes, in reverse orientation:
[kmer9756_plot.pdf](https://github.com/DorothyTamYiLing/genome_rearrangement/files/11210574/kmer9756_plot.pdf)

Inversion within genome region 1500000 and 2500000, 1500000bp boundary, kmer being intact in control genomes and split in case genomes, forward orientation kmer:
[kmer9618_plot.pdf](https://github.com/DorothyTamYiLing/genome_rearrangement/files/11210590/kmer9618_plot.pdf)

Height of arrows corresponds to proporiton of case/control genomes.

2) Plotting intact kmers without N for visualising sequence content of rearrangement :

Genome position of intact kmers without N from /ext100_merge3_ISreplaced_genomes (merging overlapping IS only) were plotted. Only kmers with unqiue genome position information (by rounding off to the nearest multiple of 1000) 

Plot of intact kmers that show rearrangements in two genome regions that are significantly associated with the structural phenotype.

![myNoNintactk_rev0fwd1](https://user-images.githubusercontent.com/34043893/231439466-2dce018c-89f4-4873-a4a2-486f655af203.png)

![myNoNintactk_rev1fwd0](https://user-images.githubusercontent.com/34043893/231439513-45c179bf-8bc6-4d1b-98b6-fe8f896ce57e.png)

Height of arrows corresponds to proporiton of case/control genomes.

(Above: 10 intact kmers that are in forward orientation in majority of structure "1" genomes, as well as reverse in majority of structure "0" genomes;
Below: 8 intact kmers that are in reverse orientation in majority of structure "1" genomes, as well as forward in majority of structure "0" genomes)

Colour indices refer to the "colour index" column in the corresponding *kmer4plot.txt file (with the same prefix), hence the corresponding kmers.

Reference: Weigand, M.R., Williams, M.M., Peng, Y., Kania, D., Pawloski, L.C., Tondella, M.L. and CDC Pertussis Working Group, 2019. Genomic survey of Bordetella pertussis diversity, United States, 2000–2013. Emerging infectious diseases, 25(4), p.780.


