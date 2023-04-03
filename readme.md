go to [tutorials](https://github.com/DorothyTamYiLing/genome_rearrangement/blob/master/tutorials.md)
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
## visualising genome rearrangements that are captured by selected significant kmers

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
 
# Pre-requisite

## Replacement of IS elements in genome set

Before using this pipeline, the repetitive elements that are speculated to have mediated rearrangement events, such as IS element, must be replaced by a short placeholder sequence (e.g. Nx15) in the genome set. This can be done using the script "iSreplace_2col.py" provided in this repository.
```
python3 iSreplace_2col.py --input <multifasta genome sequences> --coor <coordinates of IS> --out <path of output>

```
Arguments:

**genome fasta** : multifasta file of the genomes for IS replacement, gzipped, headings should be the genome IDs

**coordinates of IS** : genome coordinates of the IS element to be replaced in each genome (file format: one row per IS per genome; genome IDs in 1st column (match with the genome ID in the multifasta file), start coordinate in 2nd column, end coordinate in 3rd column; headers={sseqid	mystart	myend}, tab-delimited; example file can be found in ~/example_data/IS_coor_example.txt). Coordinates ranges of the ISs in the same gneome must not overlap. This could be the output of merge_IS.R.

**path of output** : output directory for the IS replaced genomes (one fasta per genome, not gzipped)

Example:
```
python3 iSreplace_2col.py --input ~/example_data/ISrpl_testgenomes.fasta.gz --coor ~/example_data/IS_coor_example.txt  --out path/to/your/output
```
The sequence of the IS element being replaced in this command is IS481 in _Bordetella pertussis_ genome TOHAMA1 (~/example_data/TOHAMA1_IS481_27283to28335.fasta). It is replaced in two genomes stored in ISrpl_testgenomes.fasta.gz according to the genome coordinates stored in IS_coor_example.txt.

Range of IS elements can be found in https://github.com/thanhleviet/ISfinder-sequences for multiple bacterial species.



**"Merging" IS elements (optional, recommended for gemones with high frequency of IS elements and genome rearrangements):**

It has been observed (for example in _Bordetella pertussis_) that genome rearrangements can be mediated by regions of homology that consist of more than one IS element (usually IS elements are found at the beginning and the end of the homology region). Therefore, in order to detect these type of rearangements, it is necessary to replace the whole region of homology. This can be done by merging coordinates of adjacent IS elements and treating htem as "one IS". A script "merge_IS.R" is provided for this purpose. 

```
Rscript merge_IS.R --input <coordinates of IS> --merge <integer> --extend <integer>
```
Arguments: 

**input** : genome coordinates of IS elements for merging across genomes (file format: one row per IS per genome; genome IDs in 1st column (match with the genome ID in the multifasta file), start coordinate in 2nd column, end coordinate in 3rd column; headers={sseqid	mystart	myend}, tab-delimited; example file can be found in ~/example_data/IS_coor_example.txt). Coordinates ranges of the ISs in the same gneome must not overlap.

**merge** : adjacent IS elements are merged when they are less than this number of bases apart

**extend** : for extending IS coordinates, number of bases to extend from each side of IS. This can ensure that the whole IS is replaced.

Example:
```
Rscript merge_IS.R --input ~/example_data/IS_coor_example.txt --merge 7000 --extend 100 
```

## Kmer-based GWAS

Then, a kmer-based GWAS is performed on the IS-replaced genome set (created as described above) searching for kmers that are associated with the phenotype of interested. K-mer based GWAS can be performed using pyseer (https://pyseer.readthedocs.io/en/master/index.html). K-mers can be generated using fsm-lite (https://github.com/nvalimak/fsm-lite). (See the tutorial section for detailed instructions)

From the output of pyseer, the kmers that are significantly associated with the phenotype and contain the short placeholder sequence are converted into a multi-fasta file, which is then used as one of the inputs of this pipeline (i.e. argument "kmers" of the main.sh script) for detecting potential genome rearrangement events that are associated with the phenotype of interest. (See the tutorial section for detailed instructions)
 

# Requirements:

blastn 2.6.0+, R scripting front-end version 3.4.4, Python 3.9.12

python modules: argparse, SeqIO, csv, pandas, gzip

R module: optparse
