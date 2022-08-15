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
bash main.sh <kmers> <genomes> <phenotype> <output> <size of kmers> <flanking seqeunces minimum length> <size of IS elements>

<kmers> : multifasta file of significant phenotype-associated kmers from the GWAS output (`$k_input`)

<genomes> : multifasta fie of genomes used in the GWAS (`$gen_input`)

<phenotype> : phenotype file with path (file format: sample names in 1st column, binary phenotype in 2nd column; no header, tab-delimited) (`$pheno`)

<output> : directory being created where the output files are generated (`$outdir`)

<size of kmers> : length (bp) of significant kmers (all kmers have to be the output of one GWAS hence are of the same size) (`$k_len`)

<flanking seqeunces minimum length> :  Minimum length (bp) of flanking sequences (both side) for the kmer to be blasted with the genomes; default: 30bp (`$flnk_len`)

<size of IS elements> : size of the IS element that are replaced by shorter placeholder sequence (i.e. Maximum distance (bp) between the upstream and downstream flanks in the genome for a kmer to be defined as intact kmer)  (`$flkdist`)

```
## For plotting flanks of selected kmer (visualising genome rearrangements that are captured by kmers)
```
Rscript plot_flk_kmer_prop.R --kmer <kmer ID> --phen <phenotype> --coor <myflk_behave_pheno.txt> --genome.size <genome size> --outdir <output directory> --flk.dist <flanking seqeunces minimum length>

<kmer ID> : ID of chosen kmer for plotting IS-flanking sequences

<phenotype> : phenotype file with path (file format: sample names in 1st column, binary phenotype in 2nd column; no header, tab-delimited)

<myflk_behave_pheno.txt> : myflk_behave_pheno.txt file from the output of main.sh

<genome size> : size of the genome (in thousands)

<output directory> : directory path where the plot will be generated (need to be created from before)

<flanking seqeunces minimum length> : size of the IS element that are replaced by shorter placeholder sequence (i.e. Maximum distance (bp) between the upstream and downstream flanks in the genome for a kmer to be defined as intact kmer) (same as <size of IS elements> in main.sh)

```

# Pre-requisite
Before using this pipeline, the repetitive elements that are speculated to have mediated the rearrangement events must be replaced by a short placeholder sequence (e.g. Nx15) in the genome set. Then, a kmer-based GWAS is performed searching for kmers that are associated with the phenotype of interested. The phenotype-associated kmers that contain the short placeholder sequence are picked as one of the inputs of this pipeline for detecting potential genome rearrangement events that are associated with the phenotype of interest.

# Detecting genome rearrangements

`main.sh` is the main script to run for detecting genome rearrangement.

To run `main.sh` (example):

```
bash main.sh allsig_kmer_withN.fasta 111_yearGWAS_genlist.fasta.gz  \
/home/ubuntu/Dorothy/genome_rearrangement/phenotypes.tsv \   #need to provide path for phenotype file
/home/ubuntu/Dorothy/genome_rearrangement/output 200 30 2500
```

Define all the arguments in `main.sh`:

`allsig_kmer_withN.fasta` : multifasta file of significant kmers (`$k_input`)

`111_yearGWAS_genlist.fasta` : multifasta fie of genomes (`$gen_input`)

`/home/ubuntu/Dorothy/genome_rearrangement/phenotypes.tsv` : phenotype file, no header, sample name in 1st column, binary phenotype in 2nd column, need to provide path (`$pheno`)

`/home/ubuntu/Dorothy/genome_rearrangement/output` : directory being created where the output files are generated (`$outdir`)

`200` : length (bp) of significant kmers (`$k_len`)

`30` : Minimum length (bp) of flanking sequences (both side) for the kmer to be blasted with the genomes (`$flnk_len`)

`2500` : Maximum distance (bp) between the upstream and downstream flanks in the genome for a kmer to be defined as intact kmer (`$flkdist`)


# Plotting flanks of selected kmer

`plot_flk_kmer_prop.R` is the script to run for plotting flanks of selected kmer.

To run `plot_flk_kmer_prop.R` (example):

```
Rscript plot_flk_kmer_prop.R --kmer kmer93 --phen /home/ubuntu/Dorothy/genome_rearrangement/phenotypes.tsv \
--coor /home/ubuntu/Dorothy/genome_rearrangement/output/myflk_behave_pheno.txt \
--genome.size 4000 --outdir /home/ubuntu/Dorothy/genome_rearrangement/output --flk.dist 2500
```

Define all the variables in `plot_flk_kmer_prop.R`:

`kmer`: chosen kmer for plotting flanks

`phen` : phenotype file, no header, sample name in 1st column, binary phenotype in 2nd column, need to provide path (`$pheno`)

`coor` : myflk_behave_pheno.txt file from the output

`genome.size`: size of the genome in thousands

`outdir` : path of where the plot will be generated

`flk.dist` : Maximum distance (bp) between the upstream and downstream flanks in the genome for a kmer to be defined as intact kmer (same as `$flkdist`)


# Pipeline and output files description

## Step1:

`filtering_kmer_and_blast.sh` (called by `main.sh`)

script functions:
1. filters sig. kmers for blasting by keeping only the kmers that contain flanking sequences (both side) of at least `$flnk_len` bp in size
2. blasts the filtered kmers with the genomes

output files: 
1. myout.txt (blast output file)
2. kmer_flanktooshort_4rm.txt (list of kmers that are removed due to having at least one flank being too short, i.e. <`$flnk_len`)
3. kmer_flanktooshort_flkcoor.txt (flank start and end coordinates of the kmers being removed)
4. kmer_forblast.fasta (multifasta file of kmer that have passed the filter and for blasting with genomes)



`extract_flank_coor.py` (called by `filtering_kmer_and_blast.sh`)

script functions:
1. gets the flank start and end coordinates of the sig. kmers

output files: 
1. flank_coor.txt; format: `kmerX_{end coordinate of the upstream flank}_{start coordinate of the downstream flank}_kmer size`


## Step2:

`make_flank_summary.R` (called by `main.sh`)

script functions:
1. filtering kmers based on blast hit information. Kmer passing the filter should have blast hits that fulfill the following criteria:

Criteria 1: both flanks should be found in all genomes (there should be 2 blast hits per kmer, one for each flank)

Criteria 2: Both flanks should be fully aligned with the genomes (the flank start and end coordinates in the blast hit should be consistent 
with those in flank_coor.txt)

Criteria 3: there should be no SNPs nor gaps (the values in the "mismatch" and "gap" columns should be 0 for all blast hit)

Criteria 4: Each flank should only show one unique blast hit per genome

2. for those kmers that have passed the filter, determine:

-`StartL` (genomic coordinate of the start of upstream flank)

-`EndL` (genomic coordinate of the end of upstream flank)

-`StartR` (genomic coordinate of the start of downstream flank)

-`EndR` (genomic coordinate of the end of downstream flank)

3. for those kmers that have passed the filter, determine their flank behaviours in each genome (behaviours could be 
`intact_k` (intact kmer), `mv_away` (flanks move away from each other), `swp_flk` (flanks swap in position),`mv&flp` (one flank has 
move away and flipped) according to the following rules:

To be defined as "intact kmer":

```
(StartL < EndL) & (EndL < StartR) & (StartR < EndR) & ((StartR-EndL) < flkdist)

(EndR < StartR) & (StartR < EndL) & (EndL < StartL) & ((EndL-StartR) < flkdist)
```

To be defined as "mv_away":
```

(StartL < EndL) & (EndL < StartR) & (StartR < EndR) & ((StartR-EndL) > flkdist)

(StartL > EndL) & (EndL > StartR) & (StartR > EndR) & ((EndL-StartR) > flkdist)
```

To be defined as "swp_flk":
```

(StartR < EndR) & (EndR < StartL) & (StartL < EndL)

(EndL < StartL) & (StartL < EndR) & (EndR < StartR)
```

To be defined as "mv&flp":
```

(StartL < EndL) & (EndL < EndR) & (EndR < StartR)

(EndL < StartL) & (StartL < StartR) & (StartR < EndR)

(StartR < EndR) & (EndR < EndL) & (EndL < StartL)

(EndR < StartR) & (StartR < StartL) & (StartL < EndL)
```

Flank behaviours are defined as `undefined_behave` when none of the rules above are fulfilled for that kmer.


4. for each kmer (across all genomes), count the number (and proportion) of case/control genome that in which each type of flank behaviour is found, determine the flank behaviour that is associated with case genomes and control genome respectively, include information of where in the genome  the flank behaviour take place (in form of summary statistics of genome coordinates. Finally, the most possible genome rearrangement event as indicated by each kmer is determined by the following rules:

- Define as `translocation` when case genomes/contrl genomes associated flank behaviour include `mv_away` or `swp_flk`.

- Define as `inversion` when case genomes/contrl genomes associated flank behaviour include `mv&flp`.


output files: 
1. rows_for_process.txt (blast hit of kemrs that pass the filter)
 
2. kmer_with_deletion.txt (blast hit of kmers that do not fulfill criteria 1) <sup> 1 </sup>

3. kmer_with_alignlen_issue.txt (blast hit of kmers that do not fulfill criteria 2) <sup> 1 </sup>

4. kmer_with_SNPgap.txt (blast hit of kmers that do not fulfill criteria 3) <sup> 1 </sup> 
 
5. kmer_with_multi_hits.txt (blast hit of kmers that do not fulfill criteria 4) <sup> 1 </sup>

6. myundefine_k.txt (blast hit of kmers with undefined behaviour in at least one genome) <sup> 1 </sup>

7. myflk_behave_pheno.txt (kmers with `StartL`,`EndL`,`StartR`,`EndR` and flank behaviour in each genome defined, and merged with phenotype information)

8. myall_out.txt, include the following information in columns:

**kmer**: kmer ID

**event_sum**: list of flank behaviours obserevd across genomes for this kmer, seperated by ":"

**flk_behaviour**: count and proportion of case and control genomes for each behaviour; format: `count of case genomes with behaviour/total number of case genomes (proportion): count of control genomes with behaviour/total number of control genomes (proportion)`

**case_assos**: behaviour associated with case genomes

**case_assos_prop**: proportion of case genomes with this behaviour

**ctrl_assos**: behaviour associated with control genomes

**ctrl_assos_prop**: proportion of control genomes with this behaviour

**case_assos_gp_Lflk_sumstat**: for the flank behaviours that is associated with case genomes, the summary statistics <sup> 2 </sup>  of the genome positions of the upstream flanks

**case_assos_gp_Rflk_sumstat**: for the flank behaviours that is associated with case genomes, the summary statistics <sup> 2 </sup>  of the genome positions of the downstream flanks

**ctrl_assos_gp_Lflk_sumstat**: for the flank behaviours that is associated with control genomes, the summary statistics <sup> 2 </sup> of the genome positions of the upstream flanks

**ctrl_assos_gp_Rflk_sumstat**: for the flank behaviours that is associated with control genomes, the summary statistics <sup> 2 </sup>  of the genome positions of the downstream flanks

**case_assos_gp_flkdis_sumstat**: for the flank behaviours that is associated with case genomes, the summary statistics <sup> 2 </sup>  of the distance between the upstream and downstream flanks

**ctrl_assos_gp_flkdis_sumstat**: for the flank behaviours that is associated with control genomes, the summary statistics <sup> 2 </sup>  of the distance between the upstream and downstream flanks

**event**: genome rearrangemnet event

<sup> 1 </sup> files are not produced when there is no content

<sup> 2 </sup>  summary statistics format: `minimum, 1st quantile, median, mean, 3rd quantile, maximum, standard deviation`

## Step3:
`plot_flk_kmer_prop.R`

script functions: 
For specific kmer, visualise the rearrangement event by plotting where the flanks are found in each genome
output file: 
1. kmerX_plot.pdf 
2. case_upstreamflk.txt (information <sup> 3 </sup> for plotting case upstream flank arrows in plot)
3. case_downstreamflk.txt (information <sup> 3 </sup> for plotting case downstream flank arrows in plot)
4. case_intactk.txt (information <sup> 3 </sup> for plotting case intactk arrows in plot)
5. ctrl_upstreamflk.txt (information <sup> 3 </sup> for plotting control upstream flank arrows in plot)
6. ctrl_downstreamflk.txt (information <sup> 3 </sup> for plotting control downstream flank arrows in plot)
7. ctrl_intactk.txt (information <sup> 3 </sup> for plotting control intactk arrows in plot)

<sup> 3 </sup> information includes median start and end coordinates, count and proportion of genomes
