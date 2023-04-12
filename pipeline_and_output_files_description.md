# Pipeline and output files description

## Pipeline description:

<img width="761" alt="Screenshot 2023-04-11 at 11 45 01 AM" src="https://user-images.githubusercontent.com/34043893/231137508-860801e4-6682-46e1-9f70-c2b9e98fa865.png">


## Script functions and output files

`fix_genome.py`

script functions:
1. Re-orientating whole genome assemblies so that all the assemblies start with the same chosen gene that is also in the same orientation in every genome.

output file:
1. fixed_genomes.fasta  

`merge_IS.R` (called by merge_replace_IS.sh)

script functions:

1. Merging coordinates of adjacent IS elements that are less than a certain number of bases apart (indicated by -m flag in merge_replace_IS.sh) and treating them as "one IS" during IS replacement.

2. Extend the IS cordinates for a number of base pair (indicated by -e flag) on each side to ensure complete mask of the IS. 

output files:
1. ext{}_merge{}_mergedIS.txt (new IS start and end coordinates after extension and merging)

2. ext{}_merge{}_mergedISstat.txt (statistics of IS coordinates after extension and merging, including distance between ISs and size of IS)

`iSreplace_2col.py` (called by merge_replace_IS.sh)

script functions:

1. replacing each of the "merged IS" by 15xN placeholder sequence.

output files:
1. IS-replaced genomes in a new directory (directory name: ext{}_merge{}_ISreplaced_genomes)

`class_k.py` (called by `main.sh`)

script functions:
1. separate significant kmers from GWAS into those contain N and those do not, and output them into separate fasta files for further processing. 

output files:
1. sigk_withN.fasta (kmers contain N)
2. sigk_noN.fasta (kmers do not contain N)

`filtering_kmer_and_blast.sh` (called by `main.sh`, for processing kmers containing N only)

script functions:
1. filters sig. kmers for blasting by keeping only the kmers that contain flanking sequences (both side) of at least `$flnk_len` bp in size
2. blasts the filtered kmers with the genomes

output files: 
1. myout.txt (blast output file)
2. kmer_flanktooshort_4rm.txt (list of kmers that are removed due to having at least one flank being too short, i.e. less than `$flnk_len`)
3. kmer_flanktooshort_flkcoor.txt (flank start and end coordinates of the kmers being removed)
4. kmer_forblast.fasta (multifasta file of kmer that have passed the filter and for blasting with genomes)

`extract_flank_coor.py` (called by `filtering_kmer_and_blast.sh`, for processing kmers containing N only)

script functions:
1. gets the flank start and end coordinates of the sig. kmers

output files: 
1. flank_coor.txt; format: `kmerX_{end coordinate of the upstream flank}_{start coordinate of the downstream flank}_kmer size`

`make_flank_summary.R` (called by `main.sh`, for processing kmers containing N only)

script functions:
1. filtering kmers based on blast hit information. Kmer passing the filter should have blast hits that fulfill the following criteria:

Criteria 1: the kmer should show blast match to at least 95% of the genomes

Criteria 2: both flanks of the kmers should be found in the genomes (there should be 2 blast hits per kmer, one for each flank)

Criteria 3: Each flank should only show one unique blast hit per genome

Criteria 4: Both flanks should be fully aligned with the genomes (the flank start and end coordinates in the blast hit should be consistent 
with those in flank_coor.txt)

SNPs and gaps are allowed. 

2. for those kmers that have passed the filter, determine:

-`StartL` (genomic coordinate of the start of upstream flank)

-`EndL` (genomic coordinate of the end of upstream flank)

-`StartR` (genomic coordinate of the start of downstream flank)

-`EndR` (genomic coordinate of the end of downstream flank)

3. Then, determine their flank behaviours in each genome (behaviours could be 
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

2. kmer_with_missinggenomes.txt (blast hit of kmers that do not fulfill criteria 1) <sup> 1 </sup>
 
3. kmer_with_deletion.txt (blast hit of kmers that do not fulfill criteria 2) <sup> 1 </sup>

4. kmer_with_multi_hits.txt (blast hit of kmers that do not fulfill criteria 3) <sup> 1 </sup>

5. kmer_with_alignlen_issue.txt (blast hit of kmers that do not fulfill criteria 4) <sup> 1 </sup>

6. myundefine_k.txt (blast hit of kmers with undefined behaviour in at least one genome) <sup> 1 </sup>

7. myflk_behave_pheno.txt (kmers with `StartL`,`EndL`,`StartR`,`EndR`, flank behaviour, kmer orientation (for intactk only), flank disatnce in each genome defined, and merged with phenotype information)

8. mysplitk_out.txt, include the following information in columns:

**kmer**: N-contaiing kmer ID

**event_sum**: list of flank behaviours obserevd across genomes for this kmer, seperated by ":"

**flk_behaviour**: count and proportion of case and control genomes for each behaviour; format: `count of case genomes with behaviour/total number of case genomes (proportion): count of control genomes with behaviour/total number of control genomes (proportion)`

**my0_intactk_sum**: for kmers that are intact (kmers that are not split) in control genomes, summary genome positions for the upstream and downstream flanks <sup> 2 </sup>

**my1_intactk_sum**: for kmers that are intact (kmers that are not split) in case genomes, summary genome positions for the upstream and downstream flanks <sup> 2 </sup>

**otherk**: flank behaviour other than intact k

**my0_otherk_sum**: for kmers that show flank behaviour other than intactk in control genomes, summary genome positions for the upstream and downstream flanks <sup> 2 </sup>

**my1_otherk_sum**: for kmers that show flank behaviour other than intactk in case genomes, summary genome positions for the upstream and downstream flanks <sup> 2 </sup>

**event**: genome rearrangemnet event

<sup> 1 </sup> files are not produced when there is no content

<sup> 2 </sup> summary genome positions for flanks, format: StartL_stat | StartL_sd | StartR_stat | StartR_sd | flk_dist_stat 

StartL_stat: summary statistics <sup> 3 </sup> for upstream flanks across genomes

StartL_sd: stand deviation for upstream flanks across genomes

StartR_stat: summary statistics <sup> 3 </sup> for downstream flanks across genomes

StartR_sd: standard deviation for downstream flanks across genomes

flk_dist_stat: summary statistics <sup> 3 </sup> for distance between flanks across genomes

<sup> 3 </sup>  summary statistics format: `minimum, 1st quantile, median, mean, 3rd quantile, maximum`

9. myshort_splitk_out_uniq.txt, contain kmers with unique proportion and genome position information for plotting. Include the following information in columns:

**kmer**: N-contaiing kmer ID

**intactk_mygp_ctrl_prop**: proportion of control genomes with intact kmer (round to decimal places indicated by -x flag)

**intactk_mygp_case_prop**: proportion of case genomes with intact kmer (round to decimal places indicated by -x flag)

**otherk_mygp_ctrl_prop**: proportion of control genomes with kmer of other behaviour (round to decimal places indicated by -x flag)

**otherk_mygp_case_prop**: proportion of case genomes with kmer of other behaviour (round to decimal places indicated by -x flag)

**my0_intactk_StartL_mean**: for intact kmers in control genomes, mean upstream flank start coordinate (round to decimal places indicated by -y flag)

**label**: labels generate for the kmers containing the above information 

10. myintactkwithN_out.txt, include the following information in columns:

**kmer**:	N-contaiing kmer ID

**kmer_behaviour**:	"intact k"

**flk_dist**: all observed values of flank distance across genomes

**fwdk_gen_count**:	number of genomes with kmer in forward orientation

**revk_gen_count**	: number of genomes with kmer in reverse orientation

**fwdk_0gen_prop**	: proportion of control genomes with kmer in forward orientation

**revk_0gen_prop**: proportion of control genomes with kmer in reverse orientation

**fwdk_1gen_prop**	: proportion of case genomes with kmer in forward orientation

**revk_1gen_prop**	: proportion of case genomes with kmer in reverse orientation

**fwdk_0gen_count**	: count of control genomes with kmer in forward orientation

**revk_0gen_count**: count of control genomes with kmer in reverse orientation

**fwdk_1gen_count**	: count of case genomes with kmer in forward orientation

**revk_1gen_count**	: count of case genomes with kmer in reverse orientation

**fwdk_0gen_med**	: median genome position of forward kmer in control genomes

**fwdk_0gen_sd**	: standard deviation of genome position of forward kmer in control genomes

**revk_0gen_med**: median genome position of reverser kmer in control genomes

**revk_0gen_sd**	: standard deviation of genome position of reverse kmer in control genomes

**fwdk_1gen_med**	: median genome position of forward kmer in case genomes

**fwdk_1gen_sd**	: standard deviation of genome position of forward kmer in case genomes

**revk_1gen_med**: median genome position of reverse kmer in case genomes

**revk_1gen_sd**: standard deviation of genome position of reverse kmer in case genomes


`process_sigkNoN.R` (called by main.sh)

script function:

1. blasts the intact kmers with the genomes.

2. filtering kmers based on blast hit information. Kmer passing the filter should have blast hits that fulfill the following criteria:

Criteria 1: the kmer should show blast match to at least 95% of the genomes

Criteria 2: Each flank should only show one unique blast hit per genome

Criteria 3: Both flanks should be fully aligned with the genomes (the flank start and end coordinates in the blast hit should be consistent 
with those in flank_coor.txt)

3. For intact kmers passing the fliter, make summary of their position and orientation in the genomes.

output files:

1. rows_for_process_NoN.txt (blast hit of kemrs that pass the filter)

2. kmer_with_missinggenomes_NoN.txt (blast hit of kmers that do not fulfill criteria 1) <sup> 1 </sup>
 
3. kmer_with_multi_hits_NoN.txt (blast hit of kmers that do not fulfill criteria 3) <sup> 1 </sup>

4. kmer_with_alignlen_issue_NoN.txt (blast hit of kmers that do not fulfill criteria 4) <sup> 1 </sup>

7. myflk_behave_pheno_NoN.txt (kmers with start and end genome position, kmer orientation in each genome defined, and merged with phenotype information)

8. myNoNintactk_out.txt, same column information as in myintactkwithN_out.txt

`plot_flk_kmer_prop.R` (called by main.sh)

script functions: 

For each split kmer in myshort_splitk_out_uniq.txt, visualise the rearrangement event by plotting where the flanks are found in case and control genomes

output files: (within /splitk_plots): 

1. kmerX_plot.pdf 

3. case_upstreamflk.txt (information <sup> 3 </sup> for plotting case upstream flank arrows in plot)

4. case_downstreamflk.txt (information <sup> 3 </sup> for plotting case downstream flank arrows in plot)

5. case_intactk.txt (information <sup> 3 </sup> for plotting case intactk arrows in plot)

6. ctrl_upstreamflk.txt (information <sup> 3 </sup> for plotting control upstream flank arrows in plot)

7. ctrl_downstreamflk.txt (information <sup> 3 </sup> for plotting control downstream flank arrows in plot)

8. ctrl_intactk.txt (information <sup> 3 </sup> for plotting control intactk arrows in plot)

<sup> 3 </sup> information includes median start and end coordinates, count and proportion of genomes

`plot_intactk.R` (called by main.sh)

script functions: 

1. Plotting intact kmers for visualising sequence content of rearrangement : Genome position of intact kmers were round off to the nearest 1000 and only kmers with unqiue genome position information were plotted.

output files:

1. myintactkwithN_rev1fwd0_set.txt (output of bash command before the script)
set of N-containing kmers that are in reverse orientation in majority of the case genomes and forward in orientation in majority of control genomes

2. myintactkwithN_rev0fwd1_set.txt (output of bash command before the script)
set of N-containing kmers that are in reverse orientation in majority of the control genomes and forward in orientation in majority of case genomes

3. myNoNintactk_rev1fwd0_set.txt (output of bash command before the script)
set of kmers without N that are in reverse orientation in majority of the case genomes and forward in orientation in majority of control genomes

4. myNoNintactk_rev0fwd1_set.txt (output of bash command before the script)
set of kmers without N that are in reverse orientation in majority of the control genomes and forward in orientation in majority of case genomes

5. myintactkwithN_rev0fwd1_kmer4plot.txt / myintactkwithN_rev1fwd0_kmer4plot.txt / myNoNintactk_rev0fwd1_kmer4plot.txt / myNoNintactk_rev1fwd0_kmer4plot.txt
Set of kmers with unqiue genome position <sup> 4 </sup> used in the final plot

6. myintactkwithN_rev0fwd1.png / myintactkwithN_rev1fwd0.png / myNoNintactk_rev0fwd1.png / myNoNintactk_rev1fwd0.png
Final plot of kmers set with unique genome position <sup> 4 </sup>

<sup> 4 </sup> kmers with duplicated genome positions are defined by those showing identical values after rounding off median genome position of forward kmer in control genomes (column "fwdk_0gen_med") to the closest multiplier of selected value (e.g. 100, 1000, 10000)
