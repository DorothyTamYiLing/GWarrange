#python cutgffseq_frgenome.py --gff AUSMDU00004142_NZ_CP027501.1.gff --genome AUSMDU00004142.fasta
import argparse
import pandas as pd
from Bio import SeqIO
import copy

parser = argparse.ArgumentParser()
parser.add_argument("--gff", "-f",
                    help="name of the input sig. kmer multifasta file")
parser.add_argument("--genome", "-g",
                    help="directory where the output file appear")
args = parser.parse_args()

mygff = args.gff
mygenome = args.genome


output_handle = open("homo_gffseq_frgenome.fasta", "w")

mydf = pd.read_csv(mygff,sep='\t',header=None)

mydf.columns = ["seqname","source","featuren","start","end","score","strand","frame","attribut"]

mydf. drop_duplicates(subset=["start","end"], inplace=True) 

total_rows = mydf.shape[0]

fasta_sequences = SeqIO.parse(mygenome,'fasta')
for record in fasta_sequences:
    print("read fasta")
    
for i in range(0,total_rows):
    #print(i)
    mystart=int(mydf.iloc[i,3])
    myend=int(mydf.iloc[i,4])
    if mystart > myend:
        mystart1=myend
        myend1=mystart
    else:
        mystart1=mystart
        myend1=myend
    #make a deep copy of record
    myseq=copy.deepcopy(record)
    #trim the copy only
    myseq.seq=myseq.seq[mystart1-1:myend1]
    myseq.id=str(mystart1)+"to"+str(myend1)
    myseq.description=mydf.iloc[i,8]+" len="+str(len(myseq))
    SeqIO.write(myseq, output_handle, "fasta")

