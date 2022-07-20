#Usage: python3 extract_flank_coor.py --input allsig_kmer_withN.fasta --outdir output

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

fasta_sequences = SeqIO.parse(myinput,'fasta') #this is multifasta file of the sig kmers
#fasta_sequences = SeqIO.parse("kmer_1_2.fasta",'fasta') #this is multifasta file of the sig kmers


with open(str(myoutdir)+'/flank_coor.txt', 'w') as f:
	for record in fasta_sequences:
		kmerID=record.id
		#print(kmerID)
		mykmer=record.seq
		mylen=len(mykmer)
		mystart=mykmer.index("N")  #start start coordinate of the N string
		leftflankend=mystart
		rightflankstart=mystart+15
#outputting leftflankend and rightflankstart coordinates
		f.write(str(kmerID)+"_"+str(leftflankend)+"_"+str(rightflankstart+1)+"_"+str(mylen)+'\n') 

