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

## For identify candidate repeat loci categories and estimating size of repeat sequence clusters in selected reference genome
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

**dist** : define repeat sequences to belong to the same repeat sequence cluster when they are less than this number of base pairs apart (default: 1000)


## For detecting and visualising genome rearrangement in genomes
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




