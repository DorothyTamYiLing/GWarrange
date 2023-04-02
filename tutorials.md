# Tutorial using examples input files from /example_data
 
This tutorial is based on a k-mer based GWAS conducted with pyseer using 111 American _Bordetella pertussis_ genomes as described in Weigand _et al_. 2019), with an aim of identifying genome rearrangement events that are associated with different year periods (between periods 2007-2010 and 2011-2013). 44 isolates are from year period 2007-2010 (phenotype 0) and 67 are from year period 2011-2013 (phenotype 1). Significant kmers containing the IS-elements (replaced by N x 15) from pyseer output are extracted and placed in a multifasta file (example_data/allsig_kmer_withN.fasta).

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

