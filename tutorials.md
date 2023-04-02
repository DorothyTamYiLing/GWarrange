# Tutotrial 1
 
This tutorial is based a subset of _Bordetella pertussis_ genomes as described in Weigand _et al_. 2019), in which chromosome structures were defined by exhaustive pairwise alignment. A subset of 47 genomes that display two different chromosome structures (18 genomes with structure 1 and 29 genomes with structure 2) (See Fig. 1) were selected, and a kmer-based GWAS was conducted using pyseer with an aim to identify kmers whose presence-absence patterns are associated with the chromosome structures (_i.e._ phenotype). Chromosome structure phenotype of two pairs of genomes were swapped for demonstratin purpose. 
Since it was believed that these rearrangements in genomes are likely to be mediated by homologous recombination between bl is mediated  IS481 elements in the genomes were replaced with shorter placedholder sequences (N x 15), followed 
2 Significant kmers containing the IS-elements (replaced by N x 15) from pyseer output are extracted and placed in a multifasta file (example_data/allsig_kmer_withN.fasta).

1. Detecting genome rearrangements

```
bash main.sh ~/example_data/allsig_kmer_withN.fasta ~/example_data/111_yearGWAS_genlist.fasta.gz ~/example_data/phenotypes.tsv path/to/your/output 200 30 2500
#note: full path must be provided for phenotypes.tsv file
```

Main output files (See "Pipeline and output files description" section for other outout files and more detailed description):

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

