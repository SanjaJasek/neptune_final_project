#!/usr/bin/env python

import sys
from Bio import SeqIO


# step 1: open reads file:

arguments = sys.argv
print(arguments)

ReadsFileName=sys.argv[1]

ReadsFile=open(ReadsFileName,'r')

ksize=int(sys.argv[3])

# step 2: k-merize reads

hash={}

for read in SeqIO.parse(ReadsFile, "fasta"):
	ReadSequence=str(read.seq)
	length=len(ReadSequence)
	length=length-ksize
	for num in range(0,length):
		first=num
		last=num+ksize
		kmer=ReadSequence[first:last]
		assert( len(kmer) == ksize )
		hash[kmer]=hash.get(kmer,0)+1

ReadsFile.close()

# step 3: search reference with k-mers, each match gets k-mer value added to the sequence

RefFileName=sys.argv[2]

RefFile=open(RefFileName,'r')

refseqs={}

for seq_record in SeqIO.parse(RefFile, "fasta"):

	header=seq_record.description

	seq=str(seq_record.seq)

	count = 0
	for hashkey in hash.keys():
		#print(hashkey)
		if hashkey in seq:
			#print(refseqs.get(seq))
			count = refseqs.get(header,0)+hash[hashkey]
		refseqs[header] = count

	print(str(refseqs[header])+" k-mers map to "+header)

	#refseqs[sequence][0] -> to access first number in list with key "sequence"

RefFile.close()


# TODO:
# tabular output: header, no of kmers, no of mapped reads
# mask low complexity regions
# mask is optional argument
# test different kmer sizes, how they affect mapping



