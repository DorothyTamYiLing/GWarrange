Extra requirements for the tutorials:
pyseer 1.3.10, 
unitig-caller 1.3.0, 
Frequency-based String Mining (lite)

Tip 1: To avoid files confusion, before running a new tutorial, it is advisable to remove all output files/folders from previous tutorial runs.

Tip 2: gff file of selected reference genome for generating candidate repeat loci categories stimating crepeat sequence clusters can be converted from genebank file (.gbk, downloaded from NCBI) using online tools such as http://genome2d.molgenrug.nl/g2d_tools_conversions.html


![diagrams_flowchart](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/413abac4-e692-4547-bc3b-845023cfc787)


Fig: Genome rearrangement pipeline summary chart


# Tutorial 1

This tutorial is based a subset of _Bordetella pertussis_ genomes as described in Weigand _et al_. 2019, in which chromosome structures are defined by MAUVE exhaustive pairwise alignment. A subset of 47 genomes that display two different chromosome structures (18 genomes with structure "1" and 29 genomes with structure "0" (See figure below) are used in GWAS with structure information as phenotype. Structure phenotype of two pairs of genomes are swapped for demonstration purpose.  

![tutorial1_MAUVE](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/476188bc-47e6-4390-86c8-51e60cb62146)

Fig : Two different chromosome structures that are found among 47 _Bordetella pertussis_ genomes. 

##

Go to the top level of /genome_rearrangement directory
```
cd /path/to/genome_rearrangement
```
First, selected reference genome C505 (accession: NZ_CP011687.1) is used for identifying repeat loci candidates and for estimating size of repeat loci clusters in the genome.
```
#using default parameters
bash scripts/homo_main.sh -gff ./example_data/NZ_CP011687.1_C505.gff -fna ./example_data/C505.fasta 
```
By looking at /output_homo/homo_occurence.txt, IS481 family transposase, IS481-like element IS481 family transposase and IS110-like element IS1663 family transposase are identified as most ubiqitous repeat loci categories in the reference genome. Size of largest repeat loci cluster is 5735bp (printed as standard output). 

To detect genome rearrangement associated with phenotype (-pheno), a short list of most ubiqitous repeat loci are placed in the file IS_NZ_CP025371.1.fasta (-replist) for repeat sequence detection in the input genome set (-gen). Minimal repeat sequence extension and merge parameters (-ext_mrg_min) are used for preserving rearrangement breakpoints; while a separate set of genomes are generated using extension of 7000bp (this number must be larger than the estimated size of largest repeat loci cluster, _i.e._ 5735bp) to ensure complete replacement of repeat sequence clusters by placeholder sequences, hence increasing sensitivity in detecting rearrangement boundaries (-ext_mrg_max). K-mer size of 200bp and minor allele frequency of 0.05 is used, specified as part of fsm-lite arguments (-fsmlite_arg). _gidA_ gene is used to re-oriente genomes (-startgene). Population structure is not controlled, specified as part of pyseer arguments (-pyseer_arg).

```
#concatenating genome fasta files for use
cat ./example_data/example_genomes/clus1clus2_47_genomes/*fasta.gz > ./example_data/clus1clus2_47.fna.gz

#Running GWarrange.sh. Full path should be provided to phenotype file

bash scripts/GWarrange.sh -gen ./example_data/clus1clus2_47.fna.gz -pheno /full/path/to/example_data/clus1clus2_pheno.txt \
-gen_size 4300 -startgene ./example_data/gidA.fasta -replist ./example_data/IS_NZ_CP025371.1.fasta \
-thread 8 \
-fsmlite_arg "-v -s 3 -S 44 -t tmp -m 200 -M 200" \
-pyseer_arg "--min-af 0.05 --max-af 0.95 --no-distances" \
-ext_mrg_min "100_3" -ext_mrg_max "7000_3"
```

**Visualising genome rearrangements that are captured by kmer**

From genome set with 7000bp extension and 3bp merging, split kmers are found (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes. They can be found in clus1clus2_47_ext7000_merge3_outdir/kmers_withN/mysplitk_out.txt.

1) Plotting split kmers for visualising rearrangement boundaries

Since k-mers contain highly redundant information, each k-mer is given a label, which contains its:  behaviour count and proportion in case/control genomes, genome position (represented by the mean StartL when k-mers are intact in case genomes, rounded off to two significant digits, as indicated by the -dedupk flag), and forward/reverse intact k-mers count. Only k-mers with unique labels are kept. Deduplicated k-mers can be found in output file clus1clus2_47_ext7000_merge3_outdir/kmers_withN/myshort_splitk_out_uniq.txt.

Four rearrangement boundaries are found, and they potentially refer to two inversion events, i.e. between 43000bp and 3600000bp, as well as between 1500000bp and 2500000bp, one inversion nested within the other. The four boundaries can be indicated by sixteen different significant split kmers that are mapped to each of the boundary, split in case/control genomes, and in forward/reverse orientation (plots of four split kmers are shown below as examples). Full information of these kmers can be found in output file clus1clus2_47_ext7000_merge3_outdir/kmers_withN/mysplitk_out.txt. Plots for split kmers can be found in folder clus1clus2_47_ext7000_merge3_outdir/kmers_withN/splitk_plots.

Inversion within genome region 43000 and 3600000, 43000bp boundary, kmer being intact in case genomes and split in control genomes, in forward orientation:

![kmer99_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/9a4cbda0-2637-4dea-b479-fdf1f4094f9e)


Inversion within genome region 43000 and 3600000, 3600000bp boundary, kmer being intact in control genomes and split in case genomes, in reverse orientation:

![kmer933_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/d3244de7-e602-4c4f-8ef2-a7070e4b4735)


Inversion within genome region 1500000 and 2500000, 1500000bp boundary, kmer being intact in control genomes and split in case genomes, forward orientation kmer:


![kmer932_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/5eb8f711-f103-4768-a7b7-f3050635f98c)


Inversion within genome region 1500000 and 2500000, 2500000bp boundary, kmer being intact in case genomes and split in control genomes, reverse orientation kmer:


![kmer931_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/32c29ea0-2c1b-4d22-840f-f8f6ffd51e87)


Height of arrows corresponds to proportion of case/control genomes.

2) Plotting intact kmers without N for visualising interior sequence content of rearrangements :

Genome position of intact kmers without N from /ext100_merge3_ISreplaced_genomes (minimal IS extension and merging overlapping IS only) are plotted. Only kmers with unqiue genome position information (by rounding off to the nearest multiple of 1000 using default -intkrd flag value) are kept for plotting (as shown in clus1clus2_47_ext100_merge3_outdir/kmers_noN/*kmer4plot.txt files). Plots for intact kmers can be found in folder clus1clus2_47_ext100_merge3_outdir/kmers_noN.

Plot of intact kmers that show rearrangements in two genome regions that are significantly associated with structural phenotype.

![diagrams_17](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/fe65d012-960b-457a-818c-885bc71ded60)


(Left: intact kmers that are in forward orientation in majority of structure "1" genomes, as well as in reverse orientation in majority of structure "0" genomes;
Right: intact kmers that are in reverse orientation in majority of structure "1" genomes, as well as in forward orientation in majority of structure "0" genomes)

**Important notes:**

Inversion boundaries at 1500000bp and 2500000bp are not detected using ext100_merge3_ISreplaced_genomes. It is because homologous sequence blocks/IS elements clusters at these boundaries are not completely replaced by shorter placeholder sequences. This leads to the presence of homologous sequence in flanking sequences in the kmers, hence reducing sensitivity in detecting rearrangement boundaries. To ensure complete replacement of homologous sequence blocks/IS elements clusters, we can extend the genome coordinates of each repeated sequence for a number of base pairs in both directions, and/or to merge repeated sequences that are less than a number of base pairs distance apart.

Reference: Weigand, M.R., Williams, M.M., Peng, Y., Kania, D., Pawloski, L.C., Tondella, M.L. and CDC Pertussis Working Group, 2019. Genomic survey of Bordetella pertussis diversity, United States, 2000–2013. Emerging infectious diseases, 25(4), p.780.


# Tutorial 2

_Enterococcus faecium_'s genomes are known to be enriched with IS elements, which could play important role in their genome structure's diversification (Leavis _et al._ 2007). Genome structure of 75 _Enterococcus faecium_ were characterised by socru (Page _et al._ 2020), Among which, a subset of 32 genomes displaying two different chromosome structures (21 genomes with structure "0" and 11 genomes with structure "1") (See Fig. below) are used in GWAS, with structure information as phenotype. Structure phenotype of two pairs of genomes are swapped for demonstration purpose. 

![tutorial2_MAUVE](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/7e4f08d3-ebfd-4151-a44e-c66a16f5f41a)

Fig : Two different chromosome structures found among 32 _Enterococcus faecium_ genomes. 

##

Go to the top level of /genome_rearrangement directory
```
cd /path/to/genome_rearrangement
```
First, selected reference genome AUSMDU00004142 (accession: NZ_CP027501.1) is used for identifying repeat loci candidates and for estimating size of repeat loci clusters in the genome.
```
#using default parameters
bash scripts/homo_main.sh -gff example_data/AUSMDU00004142_NZ_CP027501.1.gff -fna example_data/AUSMDU00004142.fna 
```
By looking at /output_homo/homo_occurence.txt, ISL3-like_element_ISEfa11_family_transposase, IS256-like_element_ISEf1_family_transposase and IS30_family_transposase are identified as most ubiqitous repeat loci categories in the reference genome. Size of largest repeat loci cluster is 15773bp (printed as standard output).

To detect genome rearrangement associated with phenotype (-pheno), a short list of most ubiqitous repeat loci are placed in the file Efaecium_replist.fasta (-replist) for repeat sequence detection in the input genome set (-gen). Minimal repeat sequence extension and merge parameters (-ext_mrg_min) are used for preserving rearrangement breakpoints; while a separate set of genomes are generated using extension of 17000bp (this number must be larger than the estimated size of largest repeat loci cluster, _i.e._ 15773bp) to ensure complete replacement of repeat sequence clusters by placeholder sequences, hence increasing sensitivity in detecting rearrangement boundaries (-ext_mrg_max). K-mer size of 200bp is used and only k-mers with count between 20 to 30 are output (due to the high number of k-mers produced), specified as part of fsm-lite arguments (-fsmlite_arg). _gidA_ gene is used to re-oriente genomes (-startgene). Population structure is not controlled, specified as part of pyseer arguments (-pyseer_arg). 

```
#concatenating genome fasta files for use
cat ./example_data/example_genomes/Efaecium_32genomes/*fasta.gz > ./example_data/Efaecium_32genomes.fna.gz

#Running GWarrange.sh. Full path should be provided to phenotype file

bash scripts/GWarrange.sh -gen ./example_data/Efaecium_32genomes.fna.gz -pheno /full/path/to/example_data/Efaecium32genomes_pheno_1swap.txt \
-gen_size 3000 -startgene ./example_data/dnaA.fasta -replist ./example_data/Efaecium_replist.fasta \
-thread 8 \
-fsmlite_arg "-v -t tmp -s 20 -S 30 -m 200 -M 200" \
-pyseer_arg "--print-samples --no-distances" \
-ext_mrg_min "100_3" -ext_mrg_max "17000_3" 
```
While generating the number of significant k-mers from pyseer output using 100_3 extension and merge parameters, 424,588 number of significant k-mers are found (as printed in the standard output). Due to the large number of k-mers that can be difficult to be procssed in a reasonable amount of time, it is recommended to incoporate the use of unitigs by indicating -string_type "kmers_and_unitigs". This allows unitig-callers to be run on top of fsm-lite, and the final list of significant "sequences" used for detecting genome rearrangement consisted of 1) k-mers containing placeholder seqeunces (for detecting rearrangement boundaries) and 2) unitigs without placeholder seqeunces (for detecting interior rearranged sequences).

#rerun by indicating performing both kmer-based and unitig-based GWAS
```
bash scripts/GWarrange.sh -gen ./example_data/Efaecium_32genomes.fna.gz -pheno /full/path/to/example_data/Efaecium32genomes_pheno_1swap.txt \
-gen_size 3000 -startgene ./example_data/dnaA.fasta -replist ./example_data/Efaecium_replist.fasta \
-thread 8 \
-fsmlite_arg "-v -t tmp -s 20 -S 30 -m 200 -M 200" \
-pyseer_arg "--print-samples --no-distances" \
-ext_mrg_min "100_3" -ext_mrg_max "17000_3" -string_type "kmers_and_unitigs"
```

**Visualising genome rearrangements that are captured by kmer**

Split kmers are found (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes. They can be found in Efaecium_32genomes_ext100_merge3_outdir/kmers_withN/mysplitk_out.txt.

1) Plotting split kmers for visualising rearrangement boundaries

Since kmers contain highly redundant information, only kmers with unique information (genome position, case and control count and proportion) are kept. They can be found in output file Efaecium_32genomes_ext100_merge3_outdir/kmers_withN/myshort_splitk_out_uniq.txt.

Two rearrangement boundaries are found, and they potentially refer to a single inversion event between 72000bp and 2100000bp. All of the split kmers are found to be split in case geomes and intact in control genomes. This could be explained by the absence of IS elements in the rearrangement boundaries in case genomes, which has been confirmed by manual sequence check. The two boundaries can be indicated by four different significant split kmers that are mapped to each of the boundaries in forward/reverse orientation (plots of one of these kmers are shown below). Full information of these kmers can be found in output file mysplitk_out.txt. Plots for split kmers can be found in folder Efaecium_32genomes_ext100_merge3_outdir/kmers_withN/splitk_plots.

72000bp boundary, intact kmer in forward orientation:

![kmer96750_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/6157461d-d2b1-4012-a56c-59c107247d69)


72000bp boundary, intact kmer in reverse orientation:

![kmer95738_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/f3689104-c879-4cf8-bc1c-c0d89d542486)


2100000bp boundary, intact kmer in forward orientation:

![kmer98312_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/2a94de9e-bb8c-489f-a267-5a4cfd31b5fd)


2100000bp boundary, intact kmer in reverse orientation:

![kmer97318_plot](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/dbae9891-4ec5-47ac-849f-e851481fe999)

Flanking sequences of split k-mers appeared in clusters. This is due to indels that are found in genomes leading to shifts in the rearrangement boundaries. For clarity in visualisation, arrows that are less than 40000bp apart (default value) are merged, as indicated by the --split_merge flag.

No significant split unitigs are identified as unitig-caller unitigs containing placeholder sequences are not produced by unitig-callers.

2) Plotting intact unitigs (without N) for visualising sequence content of rearrangement :

Genome positions of unitigs from /ext100_merge3_ISreplaced_genomes (minimal IS extension and merging overlapping IS only) are plotted. Only unitigs with unqiue genome position information (by rounding off to the nearest multiple of 1000, default value in -intkrd flag) are kept for plotting (as shown in Efaecium_32genomes_ext100_merge3_outdir/kmers_noN/*kmer4plot.txt files). Plots for unitigs can be found in folder Efaecium_32genomes_ext100_merge3_outdir/kmers_noN.

Plot of unitigs that show rearrangements significantly associated with structure phenotype.

![diagrams_18](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/1d5f7076-33c3-4ccb-86f8-9ff4ee8e3cdd)


(Left: intact unitigs that are in forward orientation in majority of structure "1" genomes, as well as in reverse orientation in majority of structure "0" genomes;
Right: intact unitigs that are in reverse orientation in majority of structure "1" genomes, as well as in forward orientation in majority of structure "0" genomes)


Page, A.J., Ainsworth, E.V. and Langridge, G.C., 2020. socru: typing of genome-level order and orientation around ribosomal operons in bacteria. Microbial Genomics, 6(7).

Leavis, H.L., Willems, R.J.L., van Wamel, W.J.B., Schuren, F.H., Caspers, M.P.M. and Bonten, M.J.M., 2007. Insertion sequence–driven diversification creates a globally dispersed emerging multiresistant subspecies of E. faecium. PLoS pathogens, 3(1), p.e7.

# Tutorial 3

This tutorial is based on 468 _Bordetella pertussis_ genomes with pertactin (PRN) expression information. Among them, 165 genomes show the presence of pertactin expression and 303 show absence. Pertactin expression information is taken from the supplementary material summarised in Lefrancq _et al._ 2022.

Go to the top level of /genome_rearrangement directory
```
cd /path/to/genome_rearrangement
```
Same as tutorial 1,IS481 family transposase, IS481-like element IS481 family transposase and IS110-like element IS1663 family transposase are identified as most ubiqitous repeat loci categories in the reference genome C505 (accession: NZ_CP011687.1) . Size of largest repeat loci cluster is 5735bp (printed as standard output). 
```
#using default parameters
bash script/homo_main.sh -gff example_data/C505_NZ_CP011687.1.gff -fna example_data/C505_NZ_CP011687.1.fna 
```

Population structure is controlled using phylogeny similarity matrix. Minor allele frequency of 0.05 is applied in generating k-mers and pyseer GWAS.
```
#concatenating genome fasta files for use
cat ./example_data/example_genomes/PRN_468/*fasta.gz > ./example_data/PRN_468.fna.gz

#Running GWarrange.sh. Full path should be provided to phenotype file

bash scripts/GWarrange.sh -gen ./example_data/PRN_468.fna.gz -pheno /full/path/to/example_data/prn_status_pheno.txt \
-gen_size 4300 -startgene ./example_data/gidA.fasta -replist ./example_data/IS_NZ_CP025371.1.fasta \
-thread 8 \
-pyseer_arg "--lmm --similarity /full/path/to/ClfML_kappa4.964_phylogeny_similarity.tsv --min-af 0.05 --max-af 0.95 --covariates ../example_data/covariates.txt --use-covariates 2" \
-fsmlite_arg "-v -t tmp -s 24 -S 444 -m 200 -M 200" \
-ext_mrg_min "100_3" -ext_mrg_max "7000_3"
```

**Visualising genome rearrangements that are captured by kmer**

1) Plotting split kmers for visualising rearrangement boundaries

No split kmer that indicated phenotype-associated rearrangement boundary is detected.

2) Plotting intact kmers without N for visualising interior sequence content of rearrangement :

Genome position of intact kmers without N from /ext100_merge3_ISreplaced_genomes (merging overlapping IS only) are plotted. Only kmers with unqiue genome position information (by rounding off to the nearest multiple of 1000, default value in -intkrd flag) are used in the plot. Plots for intact kmers can be found in folder PRN_468_ext100_merge3_outdir/kmers_noN.

Plots of intact kmers that show interior rearranged sequence content that are significantly associated with PRN expression phenotype.

![diagrams_19](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/b94204fa-5bf6-41b1-85b0-fb89515e314b)

(Left: intact kmers that are in forward orientation in half of PRN+ (1) genomes, as well as reverse in majority of PRN- (0) genomes;
Right: intact kmers that are in reverse orientation in half of PRN+ (1) genomes, as well as forward in majority of PRN- (0) genomes)



**Important notes:**

Some of the significant intact kmers without placeholder sequence contain sequence of pertactin autotransporter (indicated by black arrows). These kmers are not found using ext7000_merge200_ISreplaced_genomes. This is because the gene pertactin autotransporter is located immediately next to an IS element in _pertusis_ genomes, and any genome rearrangement that sits completely within the "replaced IS" region will not be detected. 

Ref: Lefrancq, N., Bouchez, V., Fernandes, N., Barkoff, A.M., Bosch, T., Dalby, T., Åkerlund, T., Darenberg, J., Fabianova, K., Vestrheim, D.F. and Fry, N.K., 2022. Global spatial dynamics and vaccine-induced fitness changes of Bordetella pertussis. Science Translational Medicine, 14(642), p.eabn3253.


# Tutorial 4

This tutorial is based on 40 simulated genomes for demonstrating flanking seqeunce behabviours during translocation. Each chromosome structure was represented by 20 genomes, and the difference in chromosome structures can be explained by two translocation events. Simulated genome sequences were taken from Bordetella pertussis genomes, and ribosomal operons that consist of mainly 16S ribosomal RNA, 23S ribosomal RNA , 5S ribosomal RNA and tRNAs were taken from Salmonella enterica genomes and inserted in each of the rearrangement boundaries produced by translocation. 

![sim_MAUVE](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/cbea1606-4af4-4390-a853-797c9e874cb4)
Fig : Schematic diagram of the two different genome structures explained by two translocation events, present in 40 simulated genomes (20 genomes with structure “0” and 20 genomes with structure “1”) 

Go to the top level of /genome_rearrangement directory
```
cd /path/to/genome_rearrangement
```
First, selected reference genome sim1.fasta is used for identifying repeat loci candidates and for estimating size of repeat loci clusters in the genome.
```
#using default parameters
bash script/homo_main.sh -gff example_data/PROKKA_sim1.gff -fna example_data/sim1.fasta 
```
By looking at /output_homo/homo_occurence.txt, IS200/IS605 family transposase, 16S ribosomal RNA, 23S ribosomal RNA and 5S ribosomal RNA  are identified as most ubiqitous repeat loci categories in the reference genome. Size of largest repeat loci cluster is 5494bp (printed as standard output).

To detect genome rearrangement associated with phenotype (-pheno), a short list of most ubiqitous repeat loci are placed in the file sim_replist.fasta (-replist) for repeat sequence detection in the input genome set (-gen). Minimal repeat sequence extension and merge parameters (-ext_mrg_min) are used for preserving rearrangement breakpoints; while a separate set of genomes are generated using extension of 7000bp (this number must be larger than the estimated size of largest repeat loci cluster, _i.e._ 5494bp) to ensure complete replacement of repeat sequence clusters by placeholder sequences, hence increasing sensitivity in detecting rearrangement boundaries (-ext_mrg_max). K-mer size of 200bp is used, specified as part of fsm-lite arguments (-fsmlite_arg). sim_startgene.fasta gene is used to re-oriente genomes (-startgene). Population structure is not controlled, specified as part of pyseer arguments (-pyseer_arg). 

```
#Concatenate genome files for use
cat ~/example_data/example_genomes/sim_translocation_39genomes/*fasta.gz > ~/example_data/sim_40genomes.fna.gz

#Running GWarrange.sh. Full path should be provided to phenotype file

bash scripts/GWarrange.sh -gen ./example_data/sim_40genomes.fna.gz -pheno /full/path/to/sim_trans_pheno.txt \
-gen_size 2532 -startgene ./example_data/sim_startgene.fasta -replist ./example_data/sim_replist.fasta \
-thread 8 -string_type "kmers" \
-fsmlite_arg "-v -t tmp -s 2 -S 38 -m 200 -M 200" \
-unitigcaller_arg "" \           
-ext_mrg_min "100_3" -ext_mrg_max "7000_3"
```

**Visualising genome rearrangements that are captured by kmer**

From genome set with 7000bp extension and 3bp merging, twelve split k-mers are found (_i.e._ flanking sequences mapped to different positions) when mapped to the original genomes. They showed mv_away and swp_flk behaviour in at least one genome. They can be found in sim_40genomes_ext7000_merge3_kmer_outdir/kmers_withN/mysplitk_out.txt.

1) Plotting split kmers for visualising rearrangement boundaries

Since k-mers contain highly redundant information, each k-mer is given a label, which contains its:  behaviour count and proportion in case/control genomes, genome position (represented by the mean StartL when k-mers are intact in case genomes, rounded off to two significant digits, as indicated by the -dedupk flag), and forward/reverse intact k-mers count. Only k-mers with unique labels are kept. Deduplicated k-mers can be found in output file sim_40genomes_ext7000_merge3_kmer_outdir/kmers_withN/myshort_splitk_out_uniq.txt.

Four rearrangement boundaries are found, and they potentially refer to two translocation events, as shown in the schematic diagram above. Full information of these kmers can be found in output file sim_40genomes_ext7000_merge3_kmer_outdir/kmers_withN/mysplitk_out.txt. Plots for split kmers can be found in folder sim_40genomes_ext7000_merge3_kmer_outdir/kmers_withN/splitk_plots.

Plots for four examples k-mers that indicate translocation boundaries is shown below:

1100K boundary, “mv_away” k-mer being intact in majority of control genomes in forward orientation and split in majority of case genomes:
![sim_kmer895](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/7707f37d-7300-4317-848e-bdc6e231fc69)

1200K boundary (in “1”/case genomes only), “swp_flk” k-mer being intact in majority of case genomes in reverse orientation and split in majority of control genomes:
![sim_kmer997](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/85ff6963-a654-4bb8-9042-ca154302a770)

1300K boundary (in “0”/control genomes only), “swp_flk” k-mer being intact in majority of control genomes in reverse orientation and split in majority of case genomes:
![sim_kmer995](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/fdfc3ac4-f424-449f-b836-3db9d360b94d)

1400K boundary, “mv_away” k-mer being intact in majority of case genomes in reverse orientation and split in majority of control genomes:
![sim_kmer99](https://github.com/DorothyTamYiLing/genome_rearrangement/assets/34043893/0158ef46-b9cc-4ffa-bc27-87074a568399)

2) Plotting intact kmers for visualising rearrangement boundaries

No significant intact k-mers without placeholder sequence was found. It is because the current method is not suitable for detecting translocations as they do not involve any change in sequence orientation.

# Tutorial 1.1 (instructions for running tutorial 1 step by step)

Go to the top level of /genome_rearrangement directory
```
cd /path/to/genome_rearrangement
```
First, selected reference genome C505 (accession: NZ_CP011687.1) is used for identifying repeat sequence categories candidates (manually stored in IS_NZ_CP025371.1.fasta) and for estimating size of repeat loci clusters in the genome (i.e. 5735bp, printed as standard output).
```
#using default parameters
bash scripts/homo_main.sh -gff ./example_data/NZ_CP011687.1_C505.gff -fna ./example_data/C505.fasta 
```

Concatenating genome fasta files for use
```
cat ./example_data/example_genomes/clus1clus2_47_genomes/*fasta.gz > ./example_data/clus1clus2_47.fna.gz
```

Unzip the genome file if neccesasry
```
gunzip example_data/clus1clus2_47.fna.gz
```

Blast starting gene gidA with genomes
```
blastn -query example_data/gidA.fasta \
-subject example_data/clus1clus2_47.fna \
-outfmt 6 -out clus1clus2_47_gidA_out.txt
```

Then, genome assemblies are re-orientated according to the position and orientation of gidA in the genomes, using the script fix_genome.py, Output file: fixed_genomes.fasta
```
python3 scripts/fix_genome.py --input example_data/clus1clus2_47.fna --mycoor clus1clus2_47_gidA_out.txt
```

Location of repeat sequence candidates in the genomes are obtained through BLAST. Multiple sequences can be placed in the same multifasta file for obtaining their genome locations in all genomes at once.
```
blastn -query example_data/IS_NZ_CP025371.1.fasta \
-subject fixed_genomes.fasta \
-outfmt 6 -out blastrep_out.txt
```

Generating repeat seqeunce-replaced genomes using specified extension and merging parameters, i.e. 7000bp in this example.
```
Rscript scripts/merge_IS.R --input blastrep_out.txt --extend 7000 --merge 3
```

Generating k-mers/unitigs using based on this set of repeat seqeunce-replaced genomes using your favourite k-mer/unitig-generation tool.

Performing GWAS for finding phenotype associated k-mers/unitigs using your favourite GWAS tools.


Build database for genome set for efficient blasting
```
makeblastdb -in clus1clus2_47.fna -dbtype nucl -out genome_db
```

Significant k-mers/unitis are then used for detecting genome rearrangement associated with phenotype
```
bash scripts/main.sh -k ext100_merge3_ISreplaced_genomes/sigk_seq.fasta \
-g example_data/clus1clus2_47.fna \
-p example_data/clus1clus2_pheno.txt -d 119360 -f 30 \
-o clus1clus2_47_ext100_merge3_outdir -s 4300 -x 2 -y 1000

bash scripts/main.sh -sigk final_sig.fasta \
-gen genome_db \
-pheno clus1clus2_pheno.txt -flk_dist 119360 -flk_len 30 \
-outdir clus1clus2_47_ext7000_merge3_outdir \
-gen_size 4300 -dedupk 2 -intkrd 1000 -thread 8 \
-exp_fac 86 -yaxis 360 -arr_dist 70 -split_h 7 -split_w 10 -merge 40000 \
-intact_h 100 -intact_w 180 -intact_res 150
```
Default values for main.sh parameters:
flk_len=30
dedupk=2
intkrd=1000
thread=8
exp_fac=86
yaxis=360
arr_dist=70
split_h=7
split_w=10
merge=40000
intact_h=100
intact_w=180
intact_res=150
