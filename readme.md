# Introduction
This pipeline detects and visualises phenotype-associated genome rearrangement events in bacterial genomes that are mediated by homologous recombination between repetitive elements (such as IS elements). 

Currently, this package only supports binary phenotypes.


# Requirements

blastn: 2.2.31+, R 4.2.2, Python 3.9.16

python modules: argparse, SeqIO, csv, pandas, gzip, os, Bio.Seq

R modules: optparse, plyr, ggplot2, ggforce, ggpubr

# Installation
genome_rearrangement can be installed using the `git clone` command
```
git clone https://github.com/DorothyTamYiLing/GWarrange.git
```

# Usage

## For identifing candidate repeat sequence categories and estimating size of repeat sequence clusters in selected reference genome
```
#example:
bash scripts/homo_main.sh -gff ./example_data/ref.gff -fna ./example_data/ref.fasta 

```
Arguments:

**gff** : gff file for selected reference genome, with fasta sequences lines removed

**fna** : complete assembly for selected reference genome

**thread_blast** : number of threads for BLAST (default: 8)

**freq** : number of occurrence/BLAST hit of a sequence in the .gff file to be defined as repeat sequence (default: 2)

**idcov** : values for -perc_identity and -qcov_hsp_perc in BLAST, in string format (default: "80_80")

**dist** : defines repeat sequences to belong to the same repeat sequence cluster when they are less than this number of base pairs apart (default: 1000)


## For detecting and visualising genome rearrangement in input genome set
```
#example:
bash scripts/GWarrange.sh -gen genomes.fna.gz -pheno /full/path/to/pheno.txt \
-gen_size 4300 -startgene startgene.fasta -replist replist.fasta \
-thread 8 \
-fsmlite_arg "-v -s 3 -S 44 -t tmp -m 200 -M 200" \
-pyseer_arg "--min-af 0.05 --max-af 0.95 --no-distances" \
-ext_mrg_min "100_3" -ext_mrg_max "7000_3"

```

Arguments:

**gen** : gzipped/gunzipped multifasta file of genomes set (original sequence without IS replacement). Must be in *.fna.gz suffix.

**pheno** : phenotype file (file format: sample names in 1st column, binary phenotype in 2nd column; no header, tab-delimited) (must be in full directory path)

**gen_size** : genome size

**startgene** : chosen gene for reorientating genomes

**replist** : representatives of repeat loci categories are aligned with reference genome

**ext_mrg_min** : minimum extending and merging neighbouring repeat sequences into blocks (default: "100_3")

**ext_mrg_max** : maximum extending and merging neighbouring repeat sequences into blocks (default: "7000_3")

**flk_len** : minimum length (bp) of flanking sequences (both side) for the kmer to be blasted with the genomes (default: 30bp)

**pyseer_arg** : additional arguments for running pyseer, apart from --phenotypes, --kmers, --output-patterns and --cpu. Full path should be provided for any additional input file. (default: "--min-af 0.05 --max-af 0.95")

**fsmlite_arg** : additional arguments for running fsm-lite, apart from -l (default: "-v -t tmp -m 200 -M 200")

**unitigcaller_arg** : additional arguments for running unitigcaller, apart from --call, --pyseer, --refs, --out and --threads. (default: "")

**string_type** : "kmer" or "kmers_and_unitigs" ("kmer" for performing k-mer-based GWAS only; "kmers_and_unitigs" for performing both k-mers and unitigs based GWAS, then instead of significant k-mers without placeholder sequences, significant unitigs will be analysed together with significant k-mer containing placeholder sequences for the purpose of efficient run time (default: "kmer")

**thread** : number of thread for BLAST, unitig-caller and pyseer (default: "8")


_Split k-mer plots parameters:_

**dedupk** : Number of significant digits (e.g. 2,3,4) for rounding off mean upstream flank start coordinate, for selecting split kmers with unique proportion and genome position information for plotting (default: 2)

**exp_fac** : how much the arrow expand horizontally for visibility in relative to the genome size. Smaller number leads to larger arrow expansion. In split kmers' plots. (default: 86)

**yaxis** : height of split kmers' plots (default: 360)

**arr_dist** : vertical distance between arrows, in split kmers' plots (default: 70)

**split_h** : height of the device in pdf function in R, in split kmers' plots (default: 7)

**split_w** : width of the device in pdf function in R, in split kmers' plots (default: 10)

**merge** : merge arrows into one when they are less than this number of base pair apart (default: 40000)

_Intact k-mer plots parameters:_

**intkrd** : The closest multiplier of selected value (e.g. 100, 1000, 10000) used for rounding off median genome position of intact kmer, for selecting intact k-mers with unique genome position information for plotting (default: 1000)

**intact_h** : height of the device in png function in R, in intact kmers' plots (default: 100)

**intact_w** : height of the device in png function in R, in intact kmers' plots (default: 180)

**intact_res** : Resolution for intact kmers' plots (default: 150)



 
For tutorials, please go to [here](https://github.com/DorothyTamYiLing/genome_rearrangement/blob/master/documents/tutorials.md) 

For pipeline and output files description, go to [here](https://github.com/DorothyTamYiLing/genome_rearrangement/blob/master/documents/pipeline_and_output_files_description.md)
 


## IS elements databaes 
Range of IS elements can be found in https://github.com/thanhleviet/ISfinder-sequences for multiple bacterial species.

## Please cite:
Tam, Y.L., Cameron, S., Preston, A. and Cowley, L., 2024. GWarrange: a pre-and post-genome-wide association studies pipeline for detecting phenotype-associated genome rearrangement events. Microbial Genomics, 10(7), p.001268.


