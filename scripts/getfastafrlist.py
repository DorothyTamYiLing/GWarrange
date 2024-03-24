import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--list", "-l",
                    help="name of the input sig. kmer multifasta file")
parser.add_argument("--out", "-o",
                    help="directory where the output file appear")
parser.add_argument("--input", "-i",
                    help="directory where the output file appear")
args = parser.parse_args()

myinput = args.input
myout = args.out
mylist = args.list

filepath=mylist
with open(filepath) as fp: 
      count=0
      sample = open(filepath,'r').read().strip().split() 
#      print(sample)

output_handle = open(myout, "w")
for i in sample:
    fasta_index = SeqIO.index(myinput, "fasta")
    record = fasta_index[str(i)]
    SeqIO.write(record, output_handle, "fasta")
