import argparse
from Bio import SeqIO
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i",
                    help="name of the core genome fasta file")
args = parser.parse_args()

myinput = args.input  #input is the fasta file being reversed

with gzip.open(str(myinput), "rt") as handle:
	fasta_sequences = SeqIO.parse(handle,'fasta')
	for record in fasta_sequences:
		print(record.id)
		with gzip.open(str(record.id)+".fasta.gz", "wt") as file:   #writing a gzip fasta file
			SeqIO.write(record, file,"fasta")
