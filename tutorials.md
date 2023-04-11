# Tutotrial 1

This tutorial is based a subset of _Bordetella pertussis_ genomes as described in Weigand _et al_. 2019), in which chromosome structures were defined by exhaustive pairwise alignment. A subset of 47 genomes that display two different chromosome structures (18 genomes with structure 1 and 29 genomes with structure 2) (See Fig. 1) were used. Structure phenotype of two pairs of genomes were swapped for demonstration purpose.  

<img width="1159" alt="Screenshot 2023-04-02 at 10 32 23 PM" src="https://user-images.githubusercontent.com/34043893/229380077-f0c15dba-7ed3-4fc3-a0d7-6d9a5f52b556.png">
Fig 1: Two different chromosome structures were found among 47 _Bordetella pertussis_ genomes. 

First, genomes asemblies from which detecting genome rearrangements are detected are re-orientated by a chosen gene, inthe case of Boredetella pertussis is gidA since it is the first gene after origin of replication. The location and orientation of gidA in the genomens are obtained by blasting it with multifasta file of genome assemblies.
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

Genome rearrangments in _Bordetella pertussis_ are belived to be largely mediated by homologous recombination of insertion sequence (IS) elements (such as IS481 and IS110). Location of IS elements () in the genomes are obtained by blasting. Sequences of more than one IS elements can be put in the same multifasta file for obtaining genome locations for all at once.
```
blastn -query ./example_data/IS_NZ_CP025371.1.fasta \
-subject fixed_genomes.fasta \
-outfmt 6 -out blastIS_out.txt
```

In addition, genome rearrangments in _Bordetella pertussis_ have also been observed to be mediated by homologous recombination of sequence blocks that consist of multiple IS elements, which are usually found at the start and end of the sequence block. These sequence blocks can be as large as several thousand base pair in size. In order to increase the chance for detecting of genome rearrangements, it is advised to replace the these sequence blocks comlpetely with placeholder sequence. 

Here, IS elements that were no more than 7000bp apart in each genome were "merged". Then, each of these "merged" IS element were replaced with shorter placedholder sequences (N x 15). Each pair of IS cordinates (start and end genome position) were extended by 5bp (default) on each side to ensure complete mask of the IS. A separaet set of IS-replaced genomes were also produced by enabling perfroming IS replacement with merging overlapping IS only (i.e. IS that are less than 3 bp apart) through passing string argument "on" to the -s flag.
```
bash ./scripts/merge_replace_IS.sh -g fixed_genomes.fasta -i blastIS_out.txt -m 7000 -s "on"
```

Each set of IS-replaced genomes using different IS merging and extending parameters were output into new directories "ext5_merge7000_ISreplaced_genomes" and "ext5_merge3_ISreplaced_genomes" (merging overlapping IS only) respectively.

Prior to GWAS, each set of IS-replaced genomes using different IS merging and extending parameters were used for generating kmers. Here, we use the kmer generation tool fsm-lite to generate kmers from genomes in directory /ext5_merge7000_ISreplaced_genomes.

```
#generating fsm-ite input file
for f in ./ext5_merge7000_ISreplaced_genomes/*_ext5_merge7000_ISreplaced.fasta; do id=$(basename "$f" _ext5_merge7000_ISreplaced.fasta); echo $id $f; done > ./ext5_merge7000_ISreplaced_genomes/clus1clus2_47_input.list

#generating kmers with size of 200 bases with minor allel frequency 0.05
fsm-lite -l ./ext5_merge7000_ISreplaced_genomes/clus1clus2_47_input.list -v -s 3 -S 44 -t tmp -m 200 -M 200 | gzip - > ./ext5_merge7000_ISreplaced_genomes/clus1clus2_47_ext5merge7000_k200_output.txt.gz
```

Then, a kmer-based GWAS was conducted using pyseer with an aim to identify kmers whose presence-absence patterns are associated with chromosome structures (_i.e._ phenotype). 

```
#adding header to phenotype file for pyseer input format
echo "samples binary" | cat - ./example_data/clus1clus2_pheno.txt > ./example_data/clus1clus2_pheno_4pyseer.txt

#run pyseer
pyseer --phenotypes ./example_data/clus1clus2_pheno_4pyseer.txt \
--kmers ./ext5_merge7000_ISreplaced_genomes/clus1clus2_47_ext5merge7000_k200_output.txt.gz \
--no-distances \
--min-af 0.05 --max-af 0.95 \
--print-samples --output-patterns ./ext5_merge7000_ISreplaced_genomes/kmer_patterns.txt \
> ./ext5_merge7000_ISreplaced_genomes/clus1clus2_47_ext5merge7000_k200_MAF0.05_nopopctrl
```

Pyseer standard output:
```
Read 47 phenotypes
Detected binary phenotype

56454 loaded variants
0 pre-filtered variants
56454 tested variants
56454 printed variants
```

Generate number of unique pattterns and p value significance threshold information:
```
/home/ubuntu/Dorothy/pyseer-master/scripts/count_patterns.py ./ext5_merge7000_ISreplaced_genomes/kmer_patterns.txt > \ ./ext5_merge7000_ISreplaced_genomes/count_pattern.txt
```
Extract kmers with p value below the the significance threshod:
```
awk '{ if ($4 <= 3.27E-04) { print } }' ./ext5_merge7000_ISreplaced_genomes/clus1clus2_47_ext5merge7000_k200_MAF0.05_nopopctrl \
> ./ext5_merge7000_ISreplaced_genomes/sigk_pyseer.txt
```

28,247 kmers were found to be significantly associated with the structural phenotype. The sequences of the kmers were extracted and placed in a multifasta file.

Extract significant kmer sequences and convert then into fasta format
```
#get the seqeunce only
awk '{print $1}' ./ext5_merge7000_ISreplaced_genomes/sigk_pyseer.txt > ./ext5_merge7000_ISreplaced_genomes/sigk_seq.txt 

number=$(cat ./ext5_merge7000_ISreplaced_genomes/sigk_seq.txt | wc -l)

rm ./ext5_merge7000_ISreplaced_genomes/header.txt   #remove existing header file, very important!

START=1
let "END=$number" 
 
for (( c=$START; c<=$END; c++ ))
do
	echo ">kmer""$c " >> ./ext5_merge7000_ISreplaced_genomes/header.txt
done

paste -d \\n ./ext5_merge7000_ISreplaced_genomes/header.txt ./ext5_merge7000_ISreplaced_genomes/sigk_seq.txt \
> ./ext5_merge7000_ISreplaced_genomes/sigk_seq.fasta
```

Then, these kmers were blasted with the original genome set for studying potential genome rearrangment that are captured by them, implemented by the following script:

```
bash ./scripts/main.sh -k ./ext5_merge7000_ISreplaced_genomes/sigk_seq.fasta \
-g ./example_data/clus1clus2_47.fna \
-p ./example_data/clus1clus2_pheno.txt -l 200 -d 70000 -f 30 \
-o ./clus1clus2_47_ext5_merge7000_tutout -s 4300

```

1009 kmers were found to be split (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes (in mysplitk_out.txt).


**Visualising genome rearrangements that are captured by kmer**

1) Plotting split kmers for visualising rearrangement boundaries

Since kmers contain highly redundant information, only kmers with unique information (genome position, case and control count and proportion) were kept. They can be found in output file "".

Four rearrangement boundaries were detected, and they potentially refer to two inversion events, one nested within the other.  Here, in total, eight kmers were detected that contained full information for the boundaries of two inversion events and they are shown below. Full information of these kmers can be found in output file "".

Inversion within genome region 43000 and 3600000, 43000bp boundary, kmer being intact in case genomes and split in control genomes, in forward orientation:

Inversion within genome region 43000 and 3600000, 43000bp boundary, reverse orientation kmer:



2) Plotting intact kmers for visualising sequence content of rearrangement :
Genome position of intact kmers were round off to the nearest 1000 and only kmers with unqiue genome position information were plotted. 

Plot of intact kmers that show rearrangements in two genome regions tha are significantly associated with the structural phenotype.

(Above: 33 intact kmers that are in forward orientation in majority of structure "1" genomes, as well as reverse in majority of structure "0" genomes;
Below: 33 intact kmers that are in reverse orientation in majority of structure "1" genomes, as well as forward in majority of structure "0" genomes)

Colour indices refer to the "colour index" column in the corresponding kmer4plot.txt file (with the same prefix), hence the corresponding kmers.

Reference: Weigand, M.R., Williams, M.M., Peng, Y., Kania, D., Pawloski, L.C., Tondella, M.L. and CDC Pertussis Working Group, 2019. Genomic survey of Bordetella pertussis diversity, United States, 2000–2013. Emerging infectious diseases, 25(4), p.780.


# Tutotrial 2

This tutorial is based on 468 _Bordetella pertussis_ genomes with pertactin (PRN) expression information. Among them, 165 genomes show the presence of pertactin expression and 303 show absence. A kmer-based GWAS was conducted using pyseer with an aim to identify kmers whose presence-absence patterns are associated with the presence of pertactiv expression. IS481s that were no more than 15000bp apart in each genome were "merged" (See "Merging IS element" in the prerequisite section in READMA.md). Then, each of these "merged" IS element were replaced with shorter placedholder sequences (N x 15), prior to the GWAS.

1,683 kmers were found to be significantly associated with the PRN expression phenotype. The sequences of the kmers were extracted and placed in a multifasta file (example_data/PRN468_merge15000_sigk.fasta).

Then, these kmers were blasted with the original genome set for studying potential genome rearrangment that are captured by them, implemented by the following script:

```
bash main.sh -k ./example_data/PRN468_merge15000_sigk.fasta \
-g ./example_data/PRN_468.fna.gz \
-p ./example_data/prn_status_pheno.txt -l 200 -d 200000 -f 30 \
-o ./example_data/PRN468_merge15000_testdir

```
1009 kmers were found to be split (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes (in mysplitk_out.txt).

**Visualising genome rearrangements that are captured by kmer**

Plotting flanks of **selected** split kmers:
```
#plotting kmer89
Rscript plot_flk_kmer_prop.R --kmer kmer89 --phen ./example_data/prn_status_pheno.txt \
--coor ./example_data/PRN468_merge15000_testdir/myflk_behave_pheno.txt \
--genome.size 4000 --outdir ./example_data/PRN468_merge15000_testdir/kmer89 --flk.dist 200000

```

<img width="969" alt="Screenshot 2023-04-02 at 11 43 22 PM" src="https://user-images.githubusercontent.com/34043893/229383098-fb929e04-7f83-4056-8b6e-a399469277c1.png">

Plotting intact kmers that DO NOT contain placeholder sequence (without N) and show significant ratio patterns:

```
#subset the kmers show show significant ratio patterns, which are reverse in majority of the 0 genomes, and forward and majority of the 1 genomes
awk -F "\t" 'NR==1; NR > 1{ if ($6 < 0.5 && $7 > 0.5 && $8 > 0.5 && $9 < 0.5) { print } }' myNoNintactk_out.txt > myNoNintactk_out_sigratio.txt

#plot the kmer set
Rscript plot_intactk.R \
--input ./example_data/PRN468_merge15000_testdir/myNoNintactk_out_sigratio.txt \
--outdir ./example_data/PRN468_merge15000_testdir \
--outname myNoNintactk_out_sigratio \
--gen_size 4300
```

<img width="822" alt="Screenshot 2023-04-02 at 11 47 02 PM" src="https://user-images.githubusercontent.com/34043893/229383210-9393c86b-138a-4aee-b4bb-0c831df907ea.png">

Colour indices refer to the "colour index" column in the corresponding kmer4plot.txt file (with the same prefix), hence the corresponding kmers.
