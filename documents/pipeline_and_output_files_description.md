# Pipeline and output files description

## Pipeline description:

<img width="553" alt="Screenshot 2023-09-18 215435" src="https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/61e1b575-5117-4a39-8124-c9d938c8e847">

## Script functions and output files

Each script has its "script functions" and "output file" section

`fix_genome.py`

script functions:
1. Re-orientating whole genome assemblies so that all the assemblies start with the same chosen gene that is also in the same orientation in every genome.

output file:
1. fixed_genomes.fasta  

###################################################

`merge_IS.R` (called by merge_replace_IS.sh)

script functions:

1. Merging coordinates of adjacent repeated sequences that are less than a certain number of bases apart (indicated by -m flag in merge_replace_IS.sh) and treating them as one "repeated sequence block" during replacement.

2. Extend the repeated sequences cordinates for a number of base pair (indicated by -e flag) on each side to ensure complete mask of repeated seqeunces. 

output files:
1. ext{}_merge{}_mergedIS.txt (new start and end coordinates after extension and merging)

2. ext{}_merge{}_mergedISstat.txt (statistics of coordinates after extension and merging, including size of and distance between repeated sequences)

###################################################

`iSreplace_2col.py` (called by merge_replace_IS.sh)

script functions:

1. replacing each of the "repeated sequence block" by 15xN placeholder sequence.

output files:
1. repeated sequence-replaced genomes in a new directory (directory name: ext{}_merge{}_ISreplaced_genomes)

###################################################

`class_k.py` (called by `main.sh`)

script functions:
1. separate significant k-mers from GWAS into those contain N and those do not, and output them into separate fasta files for further processing. 

output files:
1. sigk_withN.fasta (k-mers contain N)
2. sigk_noN.fasta (k-mers do not contain N)

###################################################

`filtering_kmer_and_blast.sh` (called by `main.sh`, for processing kmers containing N only)

script functions:
1. filters significant k-mers by keeping only the kmers that contain flanking sequences (both side) of at least `$flnk_len` bp in size
2. blasts the filtered kmers with the genomes

output files: 
1. myout.txt (blast output file)
2. kmer_flanktooshort_4rm.txt (list of kmers that are removed due to having at least one flank being too short, i.e. less than `$flnk_len`)
3. kmer_flanktooshort_flkcoor.txt (flank start and end coordinates of the kmers being removed)
4. kmer_forblast.fasta (multifasta file of kmer that have passed the filter and for blasting with genomes)

###################################################

`extract_flank_coor.py` (called by `filtering_kmer_and_blast.sh`, for processing kmers containing N only)

script functions:
1. gets the flank start and end coordinates of significant kmers

output files: 
1. flank_coor.txt; format: `kmerX_{end coordinate of the left flank}_{start coordinate of the right flank}_kmer size`

###################################################

`make_flank_summary.R` (called by `main.sh`, for processing kmers containing N only)

script functions:

1. Filtering kmers based on blast hit information. Kmer passing the filter should have blast hits that fulfill the following criteria:

Criteria 1: the k-mer should have a blast match with >=95% of samples in dataset

Criteria 2: each blast alignment should be >=90% of the flanking sequence’s length

Criteria 3: each blast alignment should show have >=95% identity match and <=10e-10 E-value

Criteria 4: the k-mer should show two blast hits to each genome, one hit for each flanking sequence

2. For those kmers that have passed the filter, determine:

-`StartL` (genome coordinate of the start of left flank)

-`EndL` (genome coordinate of the end of left flank)

-`StartR` (genome coordinate of the start of right flank)

-`EndR` (genome coordinate of the end of right flank)

3. Determine kmers' flank behaviours in each genome (behaviours could be 
`intact_k` (intact kmer), `mv_away` (flanks move away from each other), `swp_flk` (flanks swap in position),`mv&flp` (one flank has move away and flipped) according to the following rules:

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


4. For each k-mer (across all genomes), count the number (and proportion) of case/control genome that in which each type of flank behaviour is found, determine the flank behaviour that is associated with case genomes and control genome respectively, including information of where in the genome the flank behaviours take place (in the form of summary statistics of genome coordinates. Finally, the most possible genome rearrangement event as indicated by each kmer is determined by the following rules:

- Define as `translocation` when case genomes/contrl genomes associated flank behaviour include `mv_away` or `swp_flk`.

- Define as `inversion` when case genomes/contrl genomes associated flank behaviour include `mv&flp`.


output files: 
1. rows_for_process.txt (blast hit of k-mers that pass the filters) <sup> 1 </sup>

2. kmer_with_missinggenomes.txt (blast hit of k-mers that do not fulfill criteria 1) <sup> 1 </sup>

3. kmer_genomeappearonce.txt (blast hit of k-mers that do not fulfill criteria 4, k-mers showing one blast hit to at least one genome) <sup> 1 </sup>

4. kmer_with_multi_hits.txt (blast hit of k-mers that do not fulfill criteria 4, k-mers showing more than two blast hits to at least one genome) <sup> 1 </sup>

5. kmer_with_align_issue.txt (blast hit of k-mers that do not fulfill criteria 4, k-mers missing left or right hit) <sup> 1 </sup>

6. kmer_with_align_len.txt (blast hit of k-mers that do not fulfill criteria 2) <sup> 1 </sup>

7. kmer_with_ID_E_issue_k.txt (blast hit of k-mers that do not fulfill criteria 3) <sup> 1 </sup>

8. filterk_out_summary.txt (table summary fo each k-mer that do not fulfill criteria 1-4), with columns as follow:
   
column "abs_gen_k" refers to k-mers in kmer_with_missinggenomes.txt

column "del_k" refers to k-mers in kmer_genomeappearonce.txt

column "multi_hit_k" refers to k-mers in kmer_with_multi_hits.txt

column "align_issue_k" refers to k-mers in kmer_with_align_issue.txt

column "align_len_k" refers to k-mers in kmer_with_align_len.txt

column "ID_E_issue_k" refers to k-mers in kmer_with_ID_E_issue_k.txt

9. myundefine_k.txt (blast hit of k-mers with undefined behaviour in at least one genome) <sup> 1 </sup>

10. myflk_behave_pheno.txt (k-mers with `StartL`,`EndL`,`StartR`,`EndR`, flank behaviour, kmer orientation (for intactk only), flank distance in each genome defined, and merged with phenotype information)

11. mysplitk_out.txt, include the following information in columns:

**kmer**: N-contaiing k-mer ID

**event_sum**: list of flank behaviours obserevd across genomes for this k-mer, seperated by ":"

**flk_behaviour**: count and proportion of case and control genomes for each behaviour; format: `count of case genomes with behaviour/total number of case genomes (proportion): count of control genomes with behaviour/total number of control genomes (proportion)`

**my0_intactk_sum**: when kmer is intact (kmers that are not split) in control genomes, summary genome positions <sup> 2 </sup> for the left and right flanks 

**my1_intactk_sum**: when kmer is intact (kmers that are not split) in case genomes, summary genome positions <sup> 2 </sup> for the left and right flanks

**otherk**: flank behaviour other than intact k

**my0_otherk_sum**: when kmer shows flank behaviour other than intactk in control genomes, summary genome positions <sup> 2 </sup> for the left and right flanks

**my1_otherk_sum**: when kmer shows flank behaviour other than intactk in case genomes, summary genome positions <sup> 2 </sup> for the left and right flanks 

**event**: genome rearrangemnet event

<sup> 1 </sup> files are not produced when there is no content

<sup> 2 </sup> summary genome positions for flanks, format: StartL_stat | StartL_sd | StartR_stat | StartR_sd | flk_dist_stat 

StartL_stat: summary statistics <sup> 3 </sup> for left flanks across genomes

StartL_sd: stand deviation for left flanks across genomes

StartR_stat: summary statistics <sup> 3 </sup> for right flanks across genomes

StartR_sd: standard deviation for right flanks across genomes

flk_dist_stat: summary statistics <sup> 3 </sup> for distance between flanks across genomes

<sup> 3 </sup>  summary statistics format: `minimum, 1st quantile, median, mean, 3rd quantile, maximum`

12. myshort_splitk_out_uniq.txt, contain k-mers with unique behaviour count and proportion in case/control genomes, unique genome position (represented by the mean StartL when k-mers are intact in case genomes, rounded off to two significant digits, as indicated by the x flag), and unique forward/reverse intact k-mers count are used for plotting. Include the following information in columns:

**kmer**: N-contaiing k-mer ID

**intactk_mygp_ctrl_prop**: proportion of control genomes with intact k-mer 

**intactk_mygp_case_prop**: proportion of case genomes with intact k-mer 

**otherk_mygp_ctrl_prop**: proportion of control genomes with k-mer of other behaviour 

**otherk_mygp_case_prop**: proportion of case genomes with k-mer of other behaviour 

**my0_intactk_StartL_mean**: for intact k-mers in control genomes, mean left flank start coordinate, round to number of significant digits indicated by -x flag in main.sh)

**fwd_intactk_count**: count of intact k-mers in forward orientation

**rev_intactk_count**: count of intact k-mers in reverse orientation

**label**: labels generated for the k-mers containing the above information. K-mers with duplicated genome positions are defined by those showing identical values in this column

13. myintactkwithN_out.txt, include the following information in columns:

**kmer**:	N-contaiing k-mer ID

**kmer_behaviour**:	"intact k"

**flk_dist**: all observed values of flank distance across genomes

**fwdk_gen_count**:	number of genomes with k-mer in forward orientation

**revk_gen_count**	: number of genomes with k-mer in reverse orientation

**fwdk_0gen_prop**	: proportion of control genomes with k-mer in forward orientation

**revk_0gen_prop**: proportion of control genomes with k-mer in reverse orientation

**fwdk_1gen_prop**	: proportion of case genomes with k-mer in forward orientation

**revk_1gen_prop**	: proportion of case genomes with k-mer in reverse orientation

**fwdk_0gen_count**	: count of control genomes with k-mer in forward orientation

**revk_0gen_count**: count of control genomes with k-mer in reverse orientation

**fwdk_1gen_count**	: count of case genomes with k-mer in forward orientation

**revk_1gen_count**	: count of case genomes with k-mer in reverse orientation

**fwdk_0gen_med**	: median genome position of forward k-mer in control genomes

**fwdk_0gen_sd**	: standard deviation of genome position of forward k-mer in control genomes

**revk_0gen_med**: median genome position of reverser k-mer in control genomes

**revk_0gen_sd**	: standard deviation of genome position of reverse k-mer in control genomes

**fwdk_1gen_med**	: median genome position of forward k-mer in case genomes

**fwdk_1gen_sd**	: standard deviation of genome position of forward k-mer in case genomes

**revk_1gen_med**: median genome position of reverse k-mer in case genomes

**revk_1gen_sd**: standard deviation of genome position of reverse k-mer in case genomes

###################################################

`extract_knoN_length.py` (called by main.sh, for k-mers do not contain N only)

script function: 

1. extract length of k-mers that do not contain N

output file: 

1. kmernoN_length.txt

###################################################

`process_sigkNoN.R` (called by main.sh, for k-mers do not contain N only)

script function:

1. blasts the intact k-mers with the genomes.

2. filtering kmers based on blast hit information. K-mer passing the filter should have blast hits that fulfill the following criteria:

Criteria 1: the k-mer should show blast match to at least 95% of the genomes

Criteria 2: the k-mer should only show one unique blast hit per genome

Criteria 3: each blast alignment should be >=90% of the k-mer's length

Criteria 4: the k-mer should show at least 95% identity match and E-value of no more than 10e-10

3. For intact kmers passing the fliter, make summary of their position and orientation in the genomes.

output files:

1. rows_for_process_NoN.txt (blast hit of k-mers that pass the filter)

2. kmer_with_missinggenomes_NoN.txt (blast hit of k-mers that do not fulfill criteria 1) <sup> 1 </sup>
 
3. kmer_with_multi_hits_NoN.txt (blast hit of k-mers that do not fulfill criteria 2) <sup> 1 </sup>

4. kmer_with_align_len_noN.txt (blast hit of k-mers that do not fulfill criteria 3) <sup> 1 </sup>

5. kmer_with_ID_E_issue_noN.txt (blast hit of k-mers that do not fulfill criteria 4) <sup> 1 </sup>

6. filterk_out_summary_noN.txt (table summary fo each k-mer that do not fulfill criteria 1-4), with columns as follow:

column "abs_gen_k" refers to k-mers in kmer_with_missinggenomes_NoN.txt

column "multi_hit_k" refers to k-mers in kmer_with_multi_hits_NoN.txt

column "align_len_k" refers to k-mers in kmer_with_align_len_noN.txt

column "ID_E_issue_k" refers to k-mers in kmer_with_ID_E_issue_noN.txt


7. myflk_behave_pheno_NoN.txt (k-mers with start and end genome position, k-mer orientation in each genome defined, and merged with phenotype information)

8. myNoNintactk_out.txt, same column information as in myintactkwithN_out.txt

<sup> 1 </sup> files are not produced when there is no content

###################################################

`plot_flk_kmer_prop.R` (called by main.sh, for plotting split kmers)

script functions: 

1. For each split kmer in myshort_splitk_out_uniq.txt, visualise the rearrangement event by plotting where the flanks are found in case and control genomes

output files: (within output directory /splitk_plots): 

1. kmerX_plot.pdf 

3. case_upstreamflk.txt (information <sup> 4 </sup> for plotting case left flank arrows in plot)

4. case_downstreamflk.txt (information <sup> 4 </sup> for plotting case right flank arrows in plot)

5. case_intactk.txt (information <sup> 4 </sup> for plotting case intactk arrows in plot)

6. ctrl_upstreamflk.txt (information <sup> 4 </sup> for plotting control left flank arrows in plot)

7. ctrl_downstreamflk.txt (information <sup> 4 </sup> for plotting control right flank arrows in plot)

8. ctrl_intactk.txt (information <sup> 4 </sup> for plotting control intactk arrows in plot)

<sup> 4 </sup> information includes median start and end coordinates, count and proportion of genomes

###################################################

`plot_intactk.R` (called by main.sh, for plotting intact kmers)

script functions: 

1. Plotting intact kmers for visualising sequence content of rearrangement; only kmers with unqiue genome position information were plotted <sup> 5 </sup>.

output files:

1. myintactkwithN_rev1fwd0_set.txt

Set of N-containing kmers that are in reverse orientation in majority of the case genomes and forward in orientation in majority of control genomes

2. myintactkwithN_rev0fwd1_set.txt

Set of N-containing kmers that are in reverse orientation in majority of the control genomes and forward in orientation in majority of case genomes

3. myNoNintactk_rev1fwd0_set.txt

Set of kmers without N that are in reverse orientation in majority of the case genomes and forward in orientation in majority of control genomes

4. myNoNintactk_rev0fwd1_set.txt

Set of kmers without N that are in reverse orientation in majority of the control genomes and forward in orientation in majority of case genomes

5. myNoNintactk_other_set.txt

Set of kmers without N that do not belong to "rev0fwd1_set" nor "rev1fwd0_set"

6. myallintactk_rev0fwd1_set.txt

Combined set of kmers with and without N that are in reverse orientation in majority of the control genomes and forward in orientation in majority of case genomes

7. myallintactk_rev1fwd0_set.txt

Combined set of kmers with and without N that are in reverse orientation in majority of the case genomes and forward in orientation in majority of control genomes

8. myintactkwithN_rev0fwd1_kmer4plot.txt / myintactkwithN_rev1fwd0_kmer4plot.txt / myNoNintactk_rev0fwd1_kmer4plot.txt / myNoNintactk_rev1fwd0_kmer4plot.txt/ myNoNintactk_other_kmer4plot.txt / myallintactk_rev0fwd1_kmer4plot.txt / myallintactk_rev1fwd0_kmer4plot.txt

Set of kmers with unqiue genome position <sup> 5 </sup> used in the final plot

9. myintactkwithN_rev0fwd1.png / myintactkwithN_rev1fwd0.png / myNoNintactk_rev0fwd1.png / myNoNintactk_rev1fwd0.png / myNoNintactk_other.png / myallintactk_rev0fwd1.png / myallintactk_rev1fwd0.png

Final plot of kmers set with unique genome position <sup> 5 </sup>

<sup> 5 </sup> kmers with unique genome position are plotted. K-mers with duplicated genome positions are defined by those showing identical median StartL values when k-mers are in forward orientation in control genomes, after rounded off to the closest multiplier of selected value (e.g. 100, 1000, 10000), as indicated by the y flag.
