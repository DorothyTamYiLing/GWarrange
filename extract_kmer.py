#extract the kmer one by one
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--allk", "-a",
                    help="name of the core genome fasta file")                   
parser.add_argument("--kmer", "-k",
                    help="name of the core genome fasta file")
                    
args = parser.parse_args()

myinput = args.allk
myk = args.kmer

fasta_index = SeqIO.index(myinput, "fasta")
record = fasta_index[myk]
SeqIO.write(record, str(myk)+".fasta", "fasta")

