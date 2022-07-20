main.sh is the main script to run 

Step1:

filtering_kmer_and_blast.sh
script functions:
1. filters sig. kmers for blasting by keeping only the kmers that contain flanking sequences (both side) of at least $flnk_len bp in size
2. blasts the filtered kmers with the genomes
output files: 
1. myout.txt (blast output file)
2. kmer_flanktooshort_4rm.txt (list of kmers that are removed due to having at least one flank being too short, i.e. <$flnk_len)
3. kmer_flanktooshort_flkcoor.txt (flank start and end coordinates of the kmers being removed)
4. kmer_forblast.fasta (multifasta file of kmer that have passed the filter and for blasting with genomes)


extract_flank_coor.py (called by filtering_kmer_and_blast.sh)
script functions:
1. gets the flank start and end coordinates of the sig. kmers
output files: 
1. flank_coor.txt; format: kmerX_{end coordinate of the upstream flank}_{start coordinate of the upstream flank}_kmer size


Step2:

make_flank_summary.R
script functions:
1. filtering kmers based on blast hit information. Kmer passing the filter should have blast hits that fulfill the following criteria:
-both flanks should be found in all genomes (there should be 2 blast hits per kmer, one for each flank)
-Both flanks should be fully aligned with the genomes (the flank start and end coordinates in the blast hit should be consistent 
with those in flank_coor.txt)
-there should be no SNPs nor gaps (the values in the "mismatch" and "gap" columns should be 0 for all blast hit)
-Each flank should only show one unique blast hit per genome
2. for those kmers that have passed the filter, determine:
-StartL (genomic coordinate of the start of upstream flank);
-EndL (genomic coordinate of the end of upstream flank);
-StartR (genomic coordinate of the start of downstream flank);
-EndR (genomic coordinate of the end of downstream flank)
3. for those kmers that have passed the filter, determine their flank behaviours in each genome (behaviours could be 
"intact_k" (intact kmer), "mv_away" (flanks move away from each other), "swp_flk" (flanks swap in position),"mv&flp" (one flank has 
move away and flipped) or ""undefined_behave" (behaviours that are not defined according to the rules))
4. for each kmer (across all genomes), count the number (and proportion) of case/control genome that show each type of flank behaviour found, 
also include information of where in the genome the flank behaviour take place, finally give the most possible genome rearrangement event 
as indicated by each kmer.
output files: 
1. rows_for_process.txt (blast hit of kemrs that pass the filter)
2. kmer_with_deletion.txt (blast hit of kmers that do not fulfill criteria 1) *
3. kmer_with_alignlen_issue.txt (blast hit of kmers that do not fulfill criteria 2) *
4. kmer_with_SNPgap.txt (blast hit of kmers that do not fulfill criteria 3) * 
5. kmer_with_multi_hits.txt (blast hit of kmers that do not fulfill criteria 4) *
6. myundefine_k.txt (blast hit of kmers with undefined behaviour in at least one genome) *
7. myflk_behave_pheno.txt (kmers with StartL,EndL,StartR,EndR and flank behaviour in each genome defined, and merged with phenotype information)
8. myall_out.txt, include the following information in columns:
-kmer: kmer ID
-event_sum: list of flank behaviours obserevd across genomes for this kmer, seperated by ":"
-flk_behaviour: count and proportion of case and control genomes for each behaviour; format: count of case genomes with behaviour/total number of case genomes (proportion): count of control genomes with behaviour/total number of control genomes (proportion)
-case_assos: behaviour associated with case genomes
-case_assos_prop: proportion of case genomes with this behaviour
-ctrl_assos: behaviour associated with control genomes
-ctrl_assos_prop: proportion of control genomes with this behaviour
-case_assos_gp_Lflk_sumstat: for the flank behaviours that is associated with case genomes, the summary statistics** of the genome positions of the upstream flanks
-case_assos_gp_Rflk_sumstat: for the flank behaviours that is associated with case genomes, the summary statistics** of the genome positions of the downstream flanks
-ctrl_assos_gp_Lflk_sumstat: for the flank behaviours that is associated with control genomes, the summary statistics** of the genome positions of the upstream flanks
-ctrl_assos_gp_Rflk_sumstat: for the flank behaviours that is associated with control genomes, the summary statistics** of the genome positions of the downstream flanks
-case_assos_gp_flkdis_sumstat: for the flank behaviours that is associated with case genomes, the summary statistics** of the distance between the upstream and downstream flanks
-ctrl_assos_gp_flkdis_sumstat: for the flank behaviours that is associated with control genomes, the summary statistics** of the distance between the upstream and downstream flanks
-event: genome rearrangemnet event

* files are not produced when there is no content
**summary statistics format: minimum, lower quantile, mean, median, uper quantil maximum, standard deviation

Step3:
plot_flk_kmer_prop.R
script functions: 
For specific kmer, visualise the rearrangement event by plotting where the flanks are found in each genome
output file: kmerX_plot.pdf 
