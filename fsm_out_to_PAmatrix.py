#Before using this script run this in terminal to create out_kabsence.txt.gz. This specify which kmer each sample is LACKING, i.e. kmers' absence information of 
#while read line
#do
#echo $line
#x=$(zgrep -n -v ${line} fsm_out_8.txt.gz | cut -d : -f 1 | tr '\n' ' ')  #this is the gzipped fsm-lite output
#echo ${line} "," $x >> out_kabsence.txt
#done<sample_list.txt. #list of samples
#gzip out_kabsence.txt

#usage:python3 fsm_out_to_PAmatrix.py --input out_kabsence.txt.gz --kmer_count 8

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", "-i",
                    help="name genome file for IS replacement, with path")
parser.add_argument("--kmer_count", "-k",
                    help="name of file with the genomic coordinates of IS to be replaced, with path")
args = parser.parse_args()

myinput = args.input  
mykmer_count= args.kmer_count


# importing pandas
import pandas as pd

#create an empty table for storing presence absence data
#sample = open("sample_list.txt",'r').read().split('\n') 
#del sample[-1]
#for line in sample:
#    print(line) 

#load in absence data
#mytable=pd.read_table('out_kabsence.txt',delimiter=',', names = ['a', 'b'])
mytable=pd.read_csv(myinput, compression='gzip',names=['a', 'b'], sep=',')
list(mytable.columns) #get column names
print(mytable.dtypes)
#mytable['b'] = mytable['b'].map(str) #turn "b" column into string type, dont run this
sample=list(mytable['a'])
sample = [x.strip(' ') for x in sample] #remove all the white spaces

#create empty dataframe with row and column set
#mykmer_count=8 #set the number of kmer
myrow=list(range(1,int(mykmer_count)+1))
df = pd.DataFrame(columns=sample, index=myrow)
df=df.fillna(1) #fill now with all presence "1"
list(df.columns) #get column names
list(df.index) #get row names

#looping through the row to process
mytable=mytable[pd.notnull(mytable["b"])]  #remove the with NaN in 2nd column
for i in range(0,mytable.shape[0]):
#for i in range(0,4):
	myrow=i
	#print(myrow)
	#get the sample name and strip the white spaces >> this line can get all the rows with no missing
	myname=mytable.iloc[myrow,0].strip() #remove leading and ending space
	print(myname)
	#get the number and turn them into integers
	myindex=mytable.iloc[myrow,1]
	print(myindex)
	mylist=list(myindex.split(" "))
	del mylist[0]  #delete the first item in list, replace the original list
	mylist_int = list(map(int, mylist))  #convert list of string into integers
	print(mylist_int)
	#now replace with 0 for absence	
	df.loc[mylist_int, myname] = 0 #this is specifying the row index and column names
	
df.to_csv('kmer_pesence_absence.csv', sep ='\t')

#read in the file again
#new_df = pd.read_csv('kmer_pesence_absence.csv',sep="\t", index_col=0)

