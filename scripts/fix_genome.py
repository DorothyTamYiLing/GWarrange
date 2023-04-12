#Usage python3 fix_genome.py --input --coor --outdir

import argparse
import csv
import pandas as pd
from Bio import SeqIO
from Bio.Seq import reverse_complement

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i",
                    help="name genome file for IS replacement, with path")
parser.add_argument("--mycoor", "-c",
                    help="name genome file for IS replacement, with path")

args = parser.parse_args()
	
myinput = args.input  
mycoor=args.mycoor

#df = pd.read_csv("blastfirstgene_out.txt",sep='\t',header=None, names = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
df = pd.read_csv(mycoor,sep='\t',header=None, names = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
		
output_handle = open("fixed_genomes.fasta", "w")
with open(str(myinput), "rt") as handle:
	fasta_sequences = SeqIO.parse(handle,'fasta')
	for record in fasta_sequences:
		print(record.id)
		#mystart=int(df.loc[df['sseqid'] == record.id,"sstart"])
		mystart_1=df.loc[df['sseqid'] == record.id,"sstart"]
		mystart=int(mystart_1.iloc[0])
		#myend=int(df.loc[df['sseqid'] == record.id,"send"])
		myend_1=df.loc[df['sseqid'] == record.id,"send"]
		myend=int(myend_1.iloc[0])
		print("sstart : {}, send : {} ".format(int(mystart), int(myend))) 
		myseq=record.seq
		if mystart < myend and mystart>1:  #first gene in forward orientation and it is not in the beginning of assembly
			myseq=myseq[mystart-1:]+myseq[0:mystart-1]
			myseq
		elif mystart > myend:  #first gene in reverse orientation
			myseq=reverse_complement(myseq)
			myseq=myseq[(len(myseq)-mystart):]+myseq[0:(len(myseq)-mystart)]
			myseq
		record.seq=myseq
		SeqIO.write(record, output_handle, "fasta")	
