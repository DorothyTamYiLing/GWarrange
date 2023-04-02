# Tutotrial 1

This tutorial is based a subset of _Bordetella pertussis_ genomes as described in Weigand _et al_. 2019), in which chromosome structures were defined by exhaustive pairwise alignment. A subset of 47 genomes that display two different chromosome structures (18 genomes with structure 1 and 29 genomes with structure 2) (See Fig. 1) were selected, and a kmer-based GWAS was conducted using pyseer with an aim to identify kmers whose presence-absence patterns are associated with chromosome structures (_i.e._ phenotype). Structure phenotype of two pairs of genomes were swapped for demonstration purpose. 

<img width="1157" alt="Screenshot 2023-04-02 at 7 37 47 PM" src="https://user-images.githubusercontent.com/34043893/229372223-8ea48ec4-2dd4-4254-b823-581393c64b88.png">
Fig 1: Two different chromosome structures were found among 47 _Bordetella pertussis_ genomes. 

Genome rearrangments in _Bordetella pertussis_ are belived to be largely mediated by homologous recombination of IS481 elements, or of sequence blocks that consist of multiple IS elements, which are usually found at the start and end of the sequence block. These sequence blocks can be as large as several thousand bp in size. In order to increase the chance for detecting of genome rearrangements, IS481s that were no more than 7000bp apart in each genome were "merged" (See "Merging IS element" in the prerequisite section in READMA.md). Then, each of these "merged" IS element were replaced with shorter placedholder sequences (N x 15), prior to the GWAS.

14,582 kmers were found to be significantly associated with the structural phenotype. The sequences of the kmers were extracted and placed in a multifasta file (example_data/clus1clus2_sigk.fasta).

Then, these kmers were blasted with the original genome set for studying potential genome rearrangment that are captured by them, implemented by the following script:

```
bash main.sh -k ./example_data/clus1clus2_sigk_withN.fasta \
-g ./example_data/clus1clus2_47.fna.gz \
-p ./example_data/clus1clus2_pheno.txt -l 200 -d 70000 -f 30 \
-o ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir

```

Main output files (See "Pipeline and output files description" section for other outout files and more detailed description):

- sigk_withN.fasta
- sigk_noN.fasta
- myout.txt
- mynoN_out.txt
- mysplitk_out.txt
- myintactkwithN_out.txt

1009 kmers were found to be split (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes (in mysplitk_out.txt).

**Visualising genome rearrangements that are captured by kmer**

Plotting flanks of **selected** split kmers:
```
#plotting kmer9939
Rscript plot_flk_kmer_prop.R --kmer kmer9939 --phen ./example_data/clus1clus2_pheno.txt \
--coor ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir/myflk_behave_pheno.txt \
--genome.size 4000 --outdir ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir/kmer9939 --flk.dist 70000

#plotting kmer999
Rscript plot_flk_kmer_prop.R --kmer kmer999 --phen ./example_data/clus1clus2_pheno.txt \
--coor ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir/myflk_behave_pheno.txt \
--genome.size 4000 --outdir ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir/kmer999 --flk.dist 70000

```
Output file: png file contains visualisation of the rearrangement event.

<img width="882" alt="Screenshot 2023-04-02 at 9 58 16 PM" src="https://user-images.githubusercontent.com/34043893/229378714-08ace0a6-d8de-4af0-be0c-5975d2fa3170.png">

<img width="883" alt="Screenshot 2023-04-02 at 9 57 59 PM" src="https://user-images.githubusercontent.com/34043893/229378733-5af67fb3-cc1c-47c9-9f69-e6d4c294e628.png">

Plotting intact kmers:

```
#subset the kmers that are reverse in majority of the 0 genomes, and forward and majority of the 1 genomes
awk -F "\t" 'NR==1; NR > 1{ if ($6 < 0.5 && $7 > 0.5 && $8 > 0.5 && $9 < 0.5) { print } }' myintactkwithN_out.txt > myintactkwithN_rev0fwd1_set.txt

#plot the kmer set
Rscript plot_intactk.R \
--input ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir/myintactkwithN_rev0fwd1_set.txt \
--outdir ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir \
--outname myintactkwithN_rev0fwd1 \
--gen_size 4300

#subset the kmers that are reverse in majority of the 1 genomes, and forward and majority of the 0 genomes
awk -F "\t" 'NR==1; NR > 1{ if ($6 > 0.5 && $7 < 0.5 && $8 < 0.5 && $9 > 0.5) { print } }' myintactkwithN_out.txt > myintactkwithN_rev1fwd0_set.txt

#plot the kmer set
Rscript plot_intactk.R \
--input ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir/myintactkwithN_rev1fwd0_set.txt \
--outdir ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir \
--outname myintactkwithN_rev1fwd0 \

```

<img width="785" alt="Screenshot 2023-04-02 at 10 15 17 PM" src="https://user-images.githubusercontent.com/34043893/229379399-b255beac-d65f-4212-9cca-78a6447350ca.png">


<img width="720" alt="Screenshot 2023-04-02 at 10 15 05 PM" src="https://user-images.githubusercontent.com/34043893/229379384-4383bab9-b3f7-4563-8dae-0f918dcd750e.png">


Reference: Weigand, M.R., Williams, M.M., Peng, Y., Kania, D., Pawloski, L.C., Tondella, M.L. and CDC Pertussis Working Group, 2019. Genomic survey of Bordetella pertussis diversity, United States, 2000–2013. Emerging infectious diseases, 25(4), p.780.

