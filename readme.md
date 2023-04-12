# Introduction
This pipeline detects and visualises phenotype-associated genome rearrangement events in bacterial genomes that are mediated by homologous recombination between repetitive elements (such as IS elements). 

# Installation
genome_rearrangement can be installed using the `git clone` command
```
git clone https://github.com/DorothyTamYiLing/genome_rearrangement.git
```
# Usage

## For detecting genome rearrangement in genomes
```
bash main.sh -k <sigk> -g <genomes> -p <phenotype> -o <output directory> -l <size of kmers> -f <flanking seqeunces minimum length> -d <replaced size>

```

Arguments:

**k** : multifasta file of significant phenotype-associated kmers that may or may not contain repetitive elements placeholder sequence (e.g. Nx15)

**g** : gzipped multifasta file of genomes set (original sequence without IS replacement)

**p** : phenotype file (file format: sample names in 1st column, binary phenotype in 2nd column; no header, tab-delimited) 

**o** : directory being created where the output files will be generated 

**l** : length (bp) of significant kmers (all kmers have to be the output of one GWAS hence are of the same size) 

**f** :  Minimum length (bp) of flanking sequences (both side) for the kmer to be blasted with the genomes; default: 30bp 

**d** : maximum size of repetitive sequences blocks that are replaced by shorter placed holder sequence (i.e. Maximum distance (bp) between the upstream and downstream flanks in the genome for a kmer to be defined as intact kmer)  

**x** : parameter for plotting split kmers. Number of significant digits (e.g. 2,3,4; default: 2) for rounding off mean upstream flank start coordinate, for selecting split kmers with unique proportion and genome position information for plotting

**y** : parameter for plotting intact kmers. The closest multiplier of selected value (e.g. 100, 1000, 10000; default :1000) used for rounding off median genome position of intact kmer, for selecting intact kmers with unique genome position information for plotting

Example:

```
bash main.sh -k allsigk.fasta -g genomes.fna.gz -p phenotye.txt -o output_dir -l 200 -f 30 -d 200000 -x 2 -y 100 

```
 
For tutorials, please go to [here](https://github.com/DorothyTamYiLing/genome_rearrangement/blob/master/tutorials.md) 

For pipeline and output files description, go to [here](https://github.com/DorothyTamYiLing/genome_rearrangement/blob/master/pipeline%20and%20output%20files%20description.md)
 
# Pre-requisite

## Reorientating whole genome assemblies

Since most bacterial genomes are circular, genomes asemblies from which detecting genome rearrangements are detected should be re-orientated by a chosen gene, which will become the first gene in the re-orientated assemblies, in the same orientation.

First, the location and orientation of the chosen gene in the genomens are obtained by blasting it with multifasta file of genome assemblies.
```
blastn -query <fasta file of chosen gene used for genome re-orientation> \
-subject <multifasta file of genome sequences, unzipped> \
-outfmt 6 -out <output file name>

```

Then, genome assemblies are re-orientated according to the position and orientation of the chosen gene in the genomes, using the script fix_genome.py:

```
python3 ./scripts/fix_genome.py --input <multifasta genome sequences, unzipped> --mycoor <blast output file name>

```
The output file name for the genomes with same orientation is "fixed_genomes.fasta".


## Locating IS elements in the genome assemblies

Location of IS elements in the genomes are obtained by blasting. Sequences of more than one IS elements can be put in the same multifasta file for obtaining genome locations for all at once.
```
blastn -query <multifasta file for IS elements to be located in the genomes> \
-subject fixed_genomes.fasta \
-outfmt 6 -out <output file name>
```

## "Merging" IS elements Replacement of IS elements in genome set
(optional, recommended for gemones with high frequency of IS elements and genome rearrangements)

It has been observed (for example in _Bordetella pertussis_) that genome rearrangements can be mediated by regions of homology that are several thousands bp in size and consist of more than one IS element (usually IS elements are found at the beginning and the end of the homology region). Therefore, in order to detect the boundaries for these type of rearangements, it is necessary to replace the whole region of homology. This can be done by merging coordinates of adjacent IS elements (i.e. IS that are less than a certain number of bases apart) and treating them as "one IS". The merge of IS elements is optional but is recommended for gemones with high frequency of IS elements and genome rearrangements. 

For the coordinate pair (i.e. start and end position) of each IS in the genome, there is also an option to extend the IS cordinates for a number of bp on each side to ensure complete mask of the IS. Default is to extend 5bp fromm each side. 

Then, each IS element (or "merged IS") are replaced by a short placeholder sequence (e.g. Nx15) in the genome set. 

It is advised to perform IS merging with caution, as any genome rearrangement event that sits completely within the "merged IS" region will not be detected. To overcome this potential issue, user can choose to perfrom IS replacement with merging overlapping IS only (i.e. IS that are less than 3 bp apart) at the same time when users chose to merge IS that are further apart than this distance. This will lead to a separaet set of IS-replaced genomes being produced.

```
bash merge_replace_IS.sh -g fixed_genomes.fasta -i <blast outout file fo IS location in genomes> \
-e <number of bp to extend from each side of each IS, default:5> \
-m <maximum number of bp for mergeing adjacent IS, default:3 (i.e. merging overlapping IS)> \
-s <"on" or "off" string argument for outputting separate set of genomes with merging overlapping IS only>
```
Example:
```
bash ./scripts/merge_replace_IS.sh -g fixed_genomes.fasta -i IS_coor.txt -e 3 -m 5000 -s "on"
```

## Kmer-based GWAS

Then, a kmer-based GWAS is performed on the IS-replaced genome set (created as described above) searching for kmers that are associated with the phenotype of interested. K-mer based GWAS can be performed using pyseer (https://pyseer.readthedocs.io/en/master/index.html). K-mers can be generated using fsm-lite (https://github.com/nvalimak/fsm-lite). (See the tutorial section for detailed instructions)

From the output of pyseer, the kmers that are significantly associated with the phenotype and contain the short placeholder sequence are converted into a multi-fasta file, which is then used as one of the inputs of this pipeline (i.e. argument "kmers" of the main.sh script) for detecting potential genome rearrangement events that are associated with the phenotype of interest. (See the tutorial section for detailed instructions)
 
 
## IS elements databaes 
Range of IS elements can be found in https://github.com/thanhleviet/ISfinder-sequences for multiple bacterial species.



# Requirements:

blastn 2.6.0+, R scripting front-end version 3.4.4, Python 3.9.12

python modules: argparse, SeqIO, csv, pandas, gzip

R module: optparse

