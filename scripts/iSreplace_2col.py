#usage: python3 ISreplace_test.py --input E541_E537_rename.fna.gz  --coor myblastout_mergedIS.txt --out .

#install Bio in outside python
#pip3 install Bio

#script for replacing all IS elements (of different length) with 15bp "YYYYYYYYYYYYYYY" string
import argparse
import csv
import pandas as pd
from Bio import SeqIO
import gzip
import os

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i",
                    help="name genome file for IS replacement, with path")
parser.add_argument("--coor", "-c",
                    help="name of file with the genomic coordinates of IS to be replaced, with path")
parser.add_argument("--outdir", "-o",
                    help="output directory of overall scipt")

args = parser.parse_args()

myinput = args.input  
myIS= args.coor
myout= args.outdir 

mylabel=os.path.basename(str(myIS))
mylabel_1="_".join(mylabel.split("_")[0:2])

os.mkdir(str(myout))

len_rpl=15 #set the length of the replacing string
newstring="N"*len_rpl

mydf = pd.read_csv(str(myIS),sep='\t',header=0) 
#mydf = pd.read_csv("myblastout_mergedIS.txt",sep='\t',header=0)  


with open(str(myinput), "rt") as handle:
#with gzip.open(str(myinput), "rt") as handle:
#with gzip.open("E541_E537_rename.fna.gz", "rt") as handle:
	fasta_sequences = SeqIO.parse(handle,'fasta')
	for record in fasta_sequences:
		print(record.id)
#		print(repr(record.seq))
		df=mydf.loc[mydf['sseqid'] == record.id]  #select rows that refer to the sample
		total_rows = df.shape[0] 
		mystart=int(df.iloc[0,1]) #for the first IS replacement (first line coor)
#		print("mystart="+str(mystart))
		myend=int(df.iloc[0,2]) #for the first IS replacement
#		print("myend="+str(myend))
		record.seq=record.seq[:mystart-1] + newstring + record.seq[myend:] #for the first IS replacement
		lostbp_sum=0
		mylen=myend-mystart+1 #length of the IS being replaced
		lostbp=mylen-len_rpl #the numebr of bp lost after this replacement
		lostbp_sum=lostbp_sum+lostbp #update the total numebr of bp lost after replacement
		for i in range(1,total_rows):  #from the second IS replacement onwards
			mystart=int(df.iloc[i,1])
#			print("mystart="+str(mystart))
			mynewstart=mystart-lostbp_sum
#			print("mynewstart="+str(mynewstart))
			myend=int(df.iloc[i,2])
#			print("myend="+str(myend))
			mynewend=myend-lostbp_sum
#			print("mynewend="+str(mynewend))
			record.seq=record.seq[:mynewstart-1] + newstring + record.seq[mynewend:]
			mylen=myend-mystart+1
			lostbp=mylen-len_rpl #the numebr of bp lost after this replacement
			lostbp_sum=lostbp_sum+lostbp #update the total numebr of bp lost after replacement
#			print("lostbp_sum="+str(lostbp_sum))
		record.description="ISreplaced"
		SeqIO.write(record, str(myout)+"/"+str(record.id)+"_"+str(mylabel_1)+"_ISreplaced.fasta", "fasta")
		#SeqIO.write(record, str(record.id)+"_ISreplaced.fasta", "fasta")

