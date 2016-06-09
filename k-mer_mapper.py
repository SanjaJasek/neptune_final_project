#!/usr/bin/env python

import sys
from Bio import SeqIO
import math, string


arguments = sys.argv
print(arguments)

ReadsFileName=sys.argv[1]
RefFileName=sys.argv[2]
ksize=int(sys.argv[3])


def H(kseq):
	entropy = 0
	for nucl in 'ATGC':
		nuclcount=float(kseq.count(nucl))
        	p_x = nuclcount/ksize
		print('p_x for {} is {}'.format(nucl, p_x))
        	if p_x > 0:
         		entropy += - p_x*math.log(p_x, 2)
	return entropy




#################################################
# k-merize reads

ReadsFile=open(ReadsFileName,'r')

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

# calculate channon entropy for k-mers:


for hk in hash.keys():
	sh=H(hk)
	print('entropy for {} is {},'.format(hk,sh))







###############################################################
# search reference with k-mers, each match gets k-mer value added to the sequence


RefFile=open(RefFileName,'r')

refseqs={}

# printing header for output table:
print("Sequence description\tNumber of mapped k-mers\tApproximate number of mapped reads")

for seq_record in SeqIO.parse(RefFile, "fasta"):

	header=seq_record.description

	seq=str(seq_record.seq)

	count = 0
	for hashkey in hash.keys():
		#print(hashkey)
		if hashkey in seq:
			#print(refseqs.get(seq))
			count = refseqs.get(header,0)+hash[hashkey]
		refseqs[header] = float(count)

	print(header + "\t" + str(refseqs[header]) + "\t" + str(refseqs[header]/ksize) )


RefFile.close()


# TODO:
# mask low complexity regions
# mask is optional argument
# test different kmer sizes, how they affect mapping



