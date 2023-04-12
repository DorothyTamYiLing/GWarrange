# Tutotrial 1

Click into the link for downloading genome file required for tutorial:
https://drive.google.com/drive/folders/1RQU1a7kxcVSsIQr94MwSKgRf9PC0YnZO?usp=share_link

This tutorial is based a subset of _Bordetella pertussis_ genomes as described in Weigand _et al_. 2019), in which chromosome structures were defined by exhaustive pairwise alignment. A subset of 47 genomes that display two different chromosome structures (18 genomes with structure 1 and 29 genomes with structure 2) (See Fig. 1) were selected, and a kmer-based GWAS was conducted using pyseer with an aim to identify kmers whose presence-absence patterns are associated with chromosome structures (_i.e._ phenotype). Structure phenotype of two pairs of genomes were swapped for demonstration purpose. 

<img width="1159" alt="Screenshot 2023-04-02 at 10 32 23 PM" src="https://user-images.githubusercontent.com/34043893/229380077-f0c15dba-7ed3-4fc3-a0d7-6d9a5f52b556.png">
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

<img width="882" alt="Screenshot 2023-04-02 at 9 58 16 PM" src="https://user-images.githubusercontent.com/34043893/229378714-08ace0a6-d8de-4af0-be0c-5975d2fa3170.png">

<img width="883" alt="Screenshot 2023-04-02 at 9 57 59 PM" src="https://user-images.githubusercontent.com/34043893/229378733-5af67fb3-cc1c-47c9-9f69-e6d4c294e628.png">

Plotting intact kmers that contain placeholder sequence (with N):

```
#subset the intact kmers that are reverse in majority of the 0 genomes, and forward and majority of the 1 genomes
awk -F "\t" 'NR==1; NR > 1{ if ($6 < 0.5 && $7 > 0.5 && $8 > 0.5 && $9 < 0.5) { print } }' myintactkwithN_out.txt > myintactkwithN_rev0fwd1_set.txt

#plot the kmer set
Rscript plot_intactk.R \
--input ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir/myintactkwithN_rev0fwd1_set.txt \
--outdir ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir \
--outname myintactkwithN_rev0fwd1 \
--gen_size 4300

#subset the intact kmers that are reverse in majority of the 1 genomes, and forward and majority of the 0 genomes
awk -F "\t" 'NR==1; NR > 1{ if ($6 > 0.5 && $7 < 0.5 && $8 < 0.5 && $9 > 0.5) { print } }' myintactkwithN_out.txt > myintactkwithN_rev1fwd0_set.txt

#plot the kmer set
Rscript plot_intactk.R \
--input ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir/myintactkwithN_rev1fwd0_set.txt \
--outdir ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir \
--outname myintactkwithN_rev1fwd0 \

```

<img width="785" alt="Screenshot 2023-04-02 at 10 15 17 PM" src="https://user-images.githubusercontent.com/34043893/229379399-b255beac-d65f-4212-9cca-78a6447350ca.png">


<img width="720" alt="Screenshot 2023-04-02 at 10 15 05 PM" src="https://user-images.githubusercontent.com/34043893/229379384-4383bab9-b3f7-4563-8dae-0f918dcd750e.png">

Fig 2. : Plot of intact kmers that show rearrangements in two genome regions tha are significantly associated with the structural phenotype.

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
