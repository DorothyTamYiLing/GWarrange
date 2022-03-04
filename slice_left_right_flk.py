import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i",
                    help="name of the core genome fasta file")
args = parser.parse_args()

myinput = args.input
kmerID=myinput.split(".", 1)[0]
print(kmerID)

fasta_sequences = SeqIO.parse(myinput,'fasta') #this input is the output of the above script

for record in fasta_sequences:
	mykmer=record.seq
	mylen=len(mykmer)
	mystart=mykmer.index("N")  #start start coordinate of the N string
	leftflankend=mystart
	rightflankstart=mystart+15
#slicing 
	myleftflk=mykmer[0:leftflankend] 
	leftflk_len=len(myleftflk)
#slicing
	myrightflk=mykmer[rightflankstart:mylen]	
	rightflk_len=len(myrightflk)		
#outputting the flank sequences
	record.seq=myleftflk
	record.description="0:"+str(leftflankend)+":"+str(leftflk_len)
	record.id="leftflank"
	SeqIO.write(record, str(kmerID)+"_leftflank.fasta", "fasta")
#output
	record.seq=myrightflk
	record.description=str(rightflankstart+1)+":"+str(mylen)+":"+str(rightflk_len)
	record.id="rightflank"
	SeqIO.write(record, str(kmerID)+"_rightflank.fasta", "fasta")

#oututting leftflankend and rightflankstart coordinates
#with open('flank_coor.txt', 'w') as f:
#    f.write(str(leftflankend)) 
#    f.write('\n')
#    f.write(str(rightflankstart+1))
#    f.write('\n')



