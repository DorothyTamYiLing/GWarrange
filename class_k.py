#Usage: python3 class_k.py --input allsig_kmer_withN.fasta --outdir output

import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i",
                    help="name of the input sig. kmer multifasta file")
parser.add_argument("--outdir", "-o",
                    help="directory where the output file appear")
args = parser.parse_args()

myinput = args.input
myoutdir = args.outdir


mywithN = open(str(myoutdir)+"/sigk_withN.fasta", "w") 
mynoN = open(str(myoutdir)+"/sigk_noN.fasta", "w")

myrecord=SeqIO.parse(myinput, 'fasta') 
for record in myrecord:
	mycount=record.seq.count("N")
	if mycount != 0:
		SeqIO.write(record, mywithN, 'fasta')
	if mycount == 0:
		SeqIO.write(record, mynoN, 'fasta')



