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
	for num in range(0,length):
		first=num
		last=num+15
		kmer=ReadSequence[first:last]
		hash[kmer]=hash.get(kmer,0)+1

ReadsFile.close()

# step 3: search reference with k-mers, each match gets k-mer value added to the sequence
# open reference

RefFileName=sys.argv[2]

RefFile=open(RefFileName,'r')

refseqs={}

for seq_record in SeqIO.parse(RefFile, "fasta"):

	header=seq_record.id

	seq=str(seq_record.seq)
	#print(hash.keys())
	for hashkey in hash.keys():
		#print(hashkey)
		if hashkey in seq:
			#print(hash[hashkey])
			#print(refseqs.get(seq))
			refseqs[header]=refseqs.get(header,0)+hash[hashkey]
			#print(str(hash[hashkey]) +" is new hash to add")
		else:
			refseqs[header]=refseqs.get(header,0)
	print(str(refseqs[header])+" k-mers map to "+header)

	#refseqs[sequence][0] -> to access first number in list with key "sequence"




#print(refseqs)


RefFile.close()







