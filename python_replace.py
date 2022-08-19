#usage: python3 python_replace.py --start ${newstart} --end ${newend} --fasta ${fasta file to be replaced}  

#python script replacing IS with start and end coordinates
import argparse
import os
import glob
import csv
import pandas as pd
from Bio import SeqIO
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--start", "-s", 
                    help="name of the core genome fasta file")
parser.add_argument("--end","-e",
                    help="name of the core genome fasta file")
parser.add_argument("--fasta","-str",
                    help="name of the core genome fasta file")
args = parser.parse_args()

mystart = args.start
myend = args.end
myfasta = args.fasta

#print(mystart)
#print(myend)
#print(myfasta)

len_rpl=15 #set the length of the replacing string
newstring="X"*len_rpl

with open(str(myfasta), "rt") as handle:
	fasta_sequences = SeqIO.parse(handle,'fasta')
	mystrain=myfasta.split("_")[0]
	for record in fasta_sequences:
		#print(record.id)
		#print(repr(record.seq))
		record.seq=record.seq[:int(mystart)-1] + newstring + record.seq[int(myend):] #for the first IS replacement
		#forsee=record.seq[2800000:2810000]
		#print(forsee)
		record.description="ISreplaced"
		SeqIO.write(record, str(mystrain)+"_secondISrepl.fasta", "fasta")

