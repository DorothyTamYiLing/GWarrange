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
git clone https://github.com/DorothyTamYiLing/genome_rearrangement.git
```
# Usage

## For detecting genome rearrangement in genomes
```
cd /path/to/genome_rearrangement
bash scripts/main.sh -k <sigk> -g <genomes> -p <phenotype> -o <output directory> -f <flanking seqeunces minimum length> -d <replaced size>

```

Arguments:

**k** : multifasta file of significant phenotype-associated k-mers/unitigs that may or may not contain repetitive elements placeholder sequence (e.g. Nx15)

**g** : gzipped/gunzipped multifasta file of genomes set (original sequence without IS replacement)

**p** : phenotype file (file format: sample names in 1st column, binary phenotype in 2nd column; no header, tab-delimited) 

**o** : directory being created where the output files will be generated 

**f** : minimum length (bp) of flanking sequences (both side) for the kmer to be blasted with the genomes; default: 30bp 

**d** : maximum size of repetitive sequences blocks that are replaced by shorter placed holder sequence (i.e. Maximum distance (bp) between the left and right flanks in the genome for a kmer to be defined as intact kmer), written as "maxrplsize" in publication Method section  

**x** : parameter for plotting split k-mers. Number of significant digits (e.g. 2,3,4; default: 2) for rounding off mean upstream flank start coordinate, for selecting split kmers with unique proportion and genome position information for plotting

**y** : parameter for plotting intact k-mers. The closest multiplier of selected value (e.g. 100, 1000, 10000; default :1000) used for rounding off median genome position of intact kmer, for selecting intact k-mers with unique genome position information for plotting

Example:

```
cd /path/to/genome_rearrangement
bash scripts/main.sh -k allsigk.fasta -g genomes.fna.gz -p phenotye.txt -o output_dir -f 30 -d 200000 -x 2 -y 100 

```
 
For tutorials, please go to [here](https://github.com/DorothyTamYiLing/genome_rearrangement/blob/master/documents/tutorials.md) 

For pipeline and output files description, go to [here](https://github.com/DorothyTamYiLing/genome_rearrangement/blob/master/documents/pipeline_and_output_files_description.md)
 
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
python3 scripts/fix_genome.py --input <multifasta genome sequences, unzipped> --mycoor <blast output file name>

```
The output file name for the genomes with same orientation is "fixed_genomes.fasta".


## Locating repeated sequences in the genome assemblies

List of repeated sequences (such as IS elements) within the genome are replaced with short placeholder sequences of Nx15, so that flanking regions can be incorporated into the length of a k-mer. Location of repeated sequences in the genomes are obtained by blasting. Multiple repeated sequences can be put in the same multifasta file for obtaining genome locations for all at once.
```
blastn -query <multifasta file for IS elements to be located in the genomes> \
-subject fixed_genomes.fasta \
-outfmt 6 -out <output file name>
```

## "Extending and Merging" repeated seqeunces in genome set
(optional, recommended for genomes with high frequency of repeated sequences and genome rearrangements)

Repeated sequences, such as IS elements, can sometimes be found in clusters in bacterial genomes, or different types of repeated sequences can co-locate next to one another forming homologous sequence blocks. Since effective detection of genome rearrangement relies on unique mapping of flanking sequences to genomes, to ensure that flanking sequences do not contain any homologous sequence without prior information of the size of homologous sequence blocks, it would be necessary to replace the whole homologous sequence block/IS clusters by short placeholder sequences. This can be done by extending the genome coordinates of each repeated sequence for a number of base pairs in both directions, and/or to merge repeated sequences that are less than a number of base pairs distance apart .

Then, each "extended and merged" repeated sequences are replaced by a shorter placeholder sequence (e.g. Nx15) in the genome set. 

It is advised to perform extension and merging with caution, as any genome rearrangement event that sits completely within the replaced region will not be detected. To overcome this potential issue, user can choose to perfrom minimal extension (i.e. extend 100bp from each side, default) and merging overlapping repeated sequences only (i.e. those that are less than 3 bp apart, default) at the same time when users choose longer extension and/or merging repeated sequences that are further apart. This will lead to a separaet set of genomes being produced.

```
bash merge_replace_IS.sh -g fixed_genomes.fasta -i <blast outout file fo repeated seqeunces location in genomes> \
-e <number of bp to extend from each side of each IS, default:100> \
-m <maximum number of bp for mergeing adjacent repeated sequences, default:3> \
-s <"on" or "off" string argument for outputting separate set of genomes with minimal extension and merging overlapping repeated seqeunces only (i.e. using the default values)>
```
Example:
```
bash ./scripts/merge_replace_IS.sh -g fixed_genomes.fasta -i IS_coor.txt -e 7000 -m 200 -s "on"
```

## K-mer-based GWAS

Then, a kmer-based GWAS is performed on the genome set with replacement (created as described above) searching for k-mers that are associated with the phenotype of interested. K-mer-based GWAS can be performed using pyseer (https://pyseer.readthedocs.io/en/master/index.html). K-mers can be generated using fsm-lite (https://github.com/nvalimak/fsm-lite). 

From the output of pyseer, the k-mers that are significantly associated with the phenotype and contain the short placeholder sequence are converted into a multi-fasta file, which is then used as one of the inputs of this pipeline (i.e. argument "kmers" of the main.sh script) for detecting potential genome rearrangement events that are associated with the phenotype of interest. (See the tutorial example 2 for detailed instructions)

Unitigs can also be used instead of kmers, they can be generated using unitig-callers ((https://github.com/bacpop/unitig-caller) (See the tutorial section for detailed instructions)
 
## IS elements databaes 
Range of IS elements can be found in https://github.com/thanhleviet/ISfinder-sequences for multiple bacterial species.




