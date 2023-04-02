# Tutotrial 1

This tutorial is based a subset of _Bordetella pertussis_ genomes as described in Weigand _et al_. 2019), in which chromosome structures were defined by exhaustive pairwise alignment. A subset of 47 genomes that display two different chromosome structures (18 genomes with structure 1 and 29 genomes with structure 2) (See Fig. 1) were selected, and a kmer-based GWAS was conducted using pyseer with an aim to identify kmers whose presence-absence patterns are associated with chromosome structures (_i.e._ phenotype). Structure phenotype of two pairs of genomes were swapped for demonstration purpose. 

<img width="1157" alt="Screenshot 2023-04-02 at 7 37 47 PM" src="https://user-images.githubusercontent.com/34043893/229372223-8ea48ec4-2dd4-4254-b823-581393c64b88.png">
Fig 1: Two different chromosome structures were found among 47 _Bordetella pertussis_ genomes. 

Genome rearrangments in _Bordetella pertussis_ are belived to be largely mediated by homologous recombination of IS481 elements, or of sequence blocks that consist of multiple IS elements, which are usually found at the start and end of the sequence block. These sequence blocks can be as large as several thousand bp in size. In order to increase the chance for detecting of genome rearrangements, IS481s that were no more than 7000bp apart in each genome were "merged" (See "Merging IS element" in the prerequisite section in READMA.md). Then, each of these "merged" IS element were replaced with shorter placedholder sequences (N x 15), prior to the GWAS.

14,582 kmers were found to be significantly associated with the structural phenotype. The sequences of the kmers were extracted and placed in a multifasta file (example_data/clus1clus2_sigk.fasta).

Then, these kmers were blasted with the original genome set for studying potential genome rearrangment that are captured by them.

```
bash main.sh -k ./example_data/clus1clus2_sigk_withN.fasta \
-g ./example_data/clus1clus2_47.fna.gz \
-p ./example_data/clus1clus2_pheno.txt -l 200 -d 70000 -f 30 \
-o ./example_data/clus1clus2_47_merge7000GWAS_nopopctrl_testdir

```

Main output files (See "Pipeline and output files description" section for other outout files and more detailed description):

- sigk_withN.fasta
- sigk_noN.fasta
- myall_out.txt (contains key information such as rearrangement event for each kmer, genomic location of the rearrangement event, proportion of case/control genomes showing the rearrengement etc.)

2. Plotting flanks of selected kmer (visualising genome rearrangements that are captured by selected kmer)
```
Rscript plot_flk_kmer_prop.R --kmer kmer93 --phen ~/example_data/phenotypes.tsv \
--coor path/to/your/myflk_behave_pheno.txt \
--genome.size 4000 --outdir path/to/your/output --flk.dist 2500
```

Main output files (See "Pipeline and output files description" section for other outout files and more detailed description):

- A pdf file contains visualisation of the rearrangement event.


Reference: Weigand, M.R., Williams, M.M., Peng, Y., Kania, D., Pawloski, L.C., Tondella, M.L. and CDC Pertussis Working Group, 2019. Genomic survey of Bordetella pertussis diversity, United States, 2000–2013. Emerging infectious diseases, 25(4), p.780.

