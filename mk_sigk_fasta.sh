#Extracting the kmer that contain "N" from the pyseer significant kmer output file
grep "N" $1 > allsigkmer_withN.txt 

#Extract the sig kmer sequences only
awk 'NR!=1{print $1}' allsigkmer_withN.txt > sig_kmer_seq.txt

#Making the header lines for  the multi-fasta file
number=$(cat sig_kmer_seq.txt | wc -l)

START=1
let "END=$number" 
 
for (( c=$START; c<=$END; c++ ))
do
	echo ">kmer""$c " >> header.txt
done

#Create a multifasta file (input for the pipeline) of the associated kmers with N by joining the two files (header.txt and sig_kmer_list.txt) alternatively
paste -d \\n header.txt  sig_kmer_seq.txt > allsig_kmer_withN.fasta

rm header.txt
rm sig_kmer_seq.txt
rm allsigkmer_withN.txt 
