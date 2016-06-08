#!/usr/bin/env python

import sys
from Bio import SeqIO


# step 1: open reads file:

arguments = sys.argv
print(arguments)

ReadsFileName=sys.argv[1]

ReadsFile=open(ReadsFileName,'r')


# step 2: k-merize reads

hash={}

for read in SeqIO.parse(ReadsFile, "fasta"):
	ReadSequence=str(read.seq)
	length=len(ReadSequence)
	length=length-15
	print(length)
	for num in range(0,length):
		first=num
		last=num+15
		kmer=ReadSequence[first:last]
		print(kmer)
		hash[kmer]=hash.get(kmer,0)+1
print(hash)

ReadsFile.close()

# step 3: search reference with k-mers, each match gets k-mer value added to the sequence
# open reference
'''
RefFileName=sys.argv[2]

RefFile=open(RefFileName,'r')

refseqs={}
sequence=''

for seq_record in SeqIO.parse(RefFile, "fasta"):
    print(seq_record.id)
    print(seq_record.seq)




#make seq key, make header as firs element in list, number of times it appears second element


print(refseqs)


RefFile.close()


'''




