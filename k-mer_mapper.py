#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
import math, string
import matplotlib.pyplot as plt
import numpy as np
import argparse


#arguments = sys.argv
#print(arguments)
#ReadsFileName=sys.argv[1]
#RefFileName=sys.argv[2]
#ksize=int(sys.argv[3])

parser = argparse.ArgumentParser()
parser.add_argument('-reads', action="store", dest="ReadsFileName")
parser.add_argument('-ref', action="store", dest="RefFileName")
parser.add_argument('-k', action="store", dest="ksize", type=int)
#parser.add_argument('-m', action="store_true", default=False)
args = parser.parse_args()
ReadsFileName=args.ReadsFileName
RefFileName=args.RefFileName
ksize=int(args.ksize)
print(parser)



def H(kseq):
	entropy = 0
	for nucl in 'ATGC':
		nuclcount=float(kseq.count(nucl))
        	p_x = nuclcount/ksize
		#print('p_x for {} is {}'.format(nucl, p_x))
        	if p_x > 0:
         		entropy += - p_x*math.log(p_x, 2)
	return entropy

def plot(stuff):
	#print(sorted(Hdict.values()))
	plt.plot(sorted(stuff))
	plt.ylabel('Shannon entropy')
	#plt.show()
	namebase=os.path.splitext(ReadsFileName)[0]
	plt.savefig('{}_Shannon_entropy_ksize_{}.png'.format(namebase, ksize))

def kmerInfo(stuff):
	number=len(stuff)
	mean=np.mean(stuff)
	median=np.median(stuff)
	print("Number of k-mers: {}, mean Shannon entropy: {}.".format(number,mean))

def remove_lowsh(Dict):
	print("Removing low complexity k-mers.")
	for key in Dict.keys():
		val=float(Dict[key])
		#print(val)
		if val < 1:
			del Dict[key]
			del hash[key]
	
	

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

# calculate shannon entropy for k-mers:

Hdict={}
for hk in hash.keys():
	sh=H(hk)
	#print('entropy for {} is {},'.format(hk,sh))
	Hdict[hk]=sh
	#print(Hdict[hk])

plot(Hdict.values())

kmerInfo(Hdict.values())

remove_lowsh(Hdict)

plot(Hdict.values())

kmerInfo(Hdict.values())

#print(hash)

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
# has to accept fastq files for reads, not fasta
# mask is optional argument
# test different kmer sizes, how they affect mapping



