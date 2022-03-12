import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i",
                    help="name of the multifasta kmer file")
args = parser.parse_args()

myinput = args.input

fasta_sequences = SeqIO.parse(myinput,'fasta') #this is multifasta file of the sig kmers
#fasta_sequences = SeqIO.parse("kmer_1_2.fasta",'fasta') #this is multifasta file of the sig kmers


with open('flank_coor.txt', 'w') as f:
	for record in fasta_sequences:
		kmerID=record.id
		print(kmerID)
		mykmer=record.seq
		mylen=len(mykmer)
		mystart=mykmer.index("N")  #start start coordinate of the N string
		leftflankend=mystart
		rightflankstart=mystart+15
#outputting leftflankend and rightflankstart coordinates
		f.write(str(kmerID)+"_"+str(leftflankend)+"_"+str(rightflankstart+1)+"_"+str(mylen)+'\n') 



