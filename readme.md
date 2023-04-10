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

**sigk** : multifasta file of significant phenotype-associated kmers that may or may not contain repetitive elements placeholder sequence (e.g. Nx15)

**genomes** : gzipped multifasta file of genomes set (original sequence without IS replacement)

**phenotype** : phenotype file (file format: sample names in 1st column, binary phenotype in 2nd column; no header, tab-delimited) 

**output directory** : directory being created where the output files will be generated 

**size of kmers** : length (bp) of significant kmers (all kmers have to be the output of one GWAS hence are of the same size) 

**flanking seqeunces minimum length** :  Minimum length (bp) of flanking sequences (both side) for the kmer to be blasted with the genomes; default: 30bp 

**replaced size** : maximum size of repetitive sequences blocks that are replaced by shorter placed holder sequence (i.e. Maximum distance (bp) between the upstream and downstream flanks in the genome for a kmer to be defined as intact kmer)  

Example:

```
bash main.sh -k allsigk.fasta -g genomes.fna.gz -p phenotye.txt -o output_dir -l 200 -f 30 -d 200000

```
## Visualising genome rearrangements that are captured by selected significant kmers

### For plotting selected placeholder-sequence-containing kmer that is split by rearrangement (placeholder flanking sequences are plotted)
```
Rscript plot_flk_kmer_prop.R --kmer <kmer ID> --phen <phenotype> --coor <myflk_behave_pheno.txt> --genome.size <genome size> --outdir <output directory> --flk.dist <replaced size>

```

Arguments:

**kmer ID** : ID of chosen split placeholder-sequence-containing kmer for visualising genome rearrangements, for kmer IDs see first column of mysplitk_out.txt

**phenotype** : phenotype file (file format: sample names in 1st column, binary phenotype in 2nd column; no header, tab-delimited)

**myflk_behave_pheno.txt** : myflk_behave_pheno.txt file from the output of main.sh

**genome size** : size of genome (in thousands)

**output directory** : directory path where the plot will be generated

**replaced size** : maximum size of repetitive sequences blocks that are replaced by shorter placed holder sequence (i.e. Maximum distance (bp) between the upstream and downstream flanks in the genome for a kmer to be defined as intact kmer) (same as "replaced size" in main.sh)

Example:

```
Rscript plot_flk_kmer_prop.R --kmer kmer1 --phen phenotye.txt --coor myflk_behave_pheno.txt --genome.size 4000 --outdir ~/output_directory/kmer1 --flk.dist 200000
```

### For plotting selected intact kmers (containing placeholder sequence or not)
```
Rscript plot_intactk.R --input <intact k> --outdir <output directory> --outname <output name> --gen_size <genome size> 

```
Arguments:
**intact k** : selected intact kmers for visualising rearrangement. Input file containing selected rows from myintactkwithN_out.txt and/or myNoNintactk_out.txt. One line per kmer. Original header line is required.

**output directory** : directory path where the plot will be generated

**output name** : prefix of output file

**genome size** : size of genome (in thousands)

Example:

```
Rscript plot_intactk.R --input myNoNintactk_out_selected.txt --outdir ~/output_dir \
--outname myNoNintactk_out_selected \
--gen_size 4300

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
python3 fix_genome.py --input <multifasta genome sequences, unzipped> --mycoor <blast output file name> --outdir <output directory>

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

It has been observed (for example in _Bordetella pertussis_) that genome rearrangements can be mediated by regions of homology that consist of more than one IS element (usually IS elements are found at the beginning and the end of the homology region). Therefore, in order to detect these type of rearangements, it is necessary to replace the whole region of homology. This can be done by merging coordinates of adjacent IS elements and treating htem as "one IS". The merge of IS elements is optional but is recommended for gemones with high frequency of IS elements and genome rearrangements. Then, the IS/mergedIS are replaced by a short placeholder sequence (e.g. Nx15) in the genome set

```
bash merge_replace_IS.sh -g fixed_genomes.fasta -i <blast outout file fo IS location in genomes> \
-e <number of bp to extend from each side of each IS> -m <maximum number of bp for mergeing adjacent IS> \
-s <"on" or "off" string argument for outputting genomes with no extension and merging overlapping IS only>
```
Example:
```
bash merge_replace_IS.sh -g fixed_genomes.fasta -i IS_coor.txt -e 10 -m 5000 -s "on"
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
