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

RefFileName=sys.argv[2]

RefFile=open(RefFileName,'r')

refseqs={}

for seq_record in SeqIO.parse(RefFile, "fasta"):
	#print(seq_record.id)
	header=seq_record.id
	#print(seq_record.seq)
	seq=str(seq_record.seq)
	#print(hash.keys())
	for hashkey in hash.keys():
		print(hashkey)
		if hashkey in seq:
			#print(hash[hashkey])
			#print('printing')
			#print(refseqs.get(seq))
			adding=0
			if seq in refseqs:
				adding=refseqs[seq][1]
				print(str(adding)+" alreaddy mapped")
			refseqs[seq]=[header, adding + hash[hashkey]]
			print(str(hash[hashkey]) +" is new hash to add")
			print(refseqs[seq])

	#refseqs[sequence][0] -> to access readnumber




#make seq key, make header as first element in list, number of times it appears second element


#print(refseqs)


RefFile.close()







