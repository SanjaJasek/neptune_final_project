#!/usr/bin/env python

import sys, os, math, argparse
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np


#arguments = sys.argv
#print(arguments)
#ReadsFileName=sys.argv[1]
#RefFileName=sys.argv[2]
#ksize=int(sys.argv[3])

parser = argparse.ArgumentParser()
parser.add_argument('-reads', action="store", dest="ReadsFileName", help="file with reads")
parser.add_argument('-ref', action="store", dest="RefFileName", help="file with reference")
parser.add_argument('-k', action="store", dest="ksize", type=int, help="k-mer size")
parser.add_argument('-m', action="store_true", default=False, help="mask low complexity regions")
#args = parser.parse_args()
#ReadsFileName=args.ReadsFileName
#RefFileName=args.RefFileName
#ksize=int(args.ksize)
#print(parser)

try:
	args = parser.parse_args()
	ReadsFileName=args.ReadsFileName
	RefFileName=args.RefFileName
	ksize=int(args.ksize)
except:
	parser.print_help()
	sys.exit(0)


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
	plt.savefig('{}_Shannon_entropy_ksize_{}.png'.format(namebase, ksize))

def kmerInfo(stuff):
	number=len(stuff)
	mean=np.mean(stuff)
	median=np.median(stuff)
	statout.write("\nSTATISTICS\n\nNumber of k-mers: {}\nMean Shannon entropy: {}.".format(number,mean))

def remove_lowsh(Dict):
	#print("Removing low complexity k-mers.")
	for key in Dict.keys():
		val=float(Dict[key])
		#print(val)
		if val < 1:
			del Dict[key]
			del hash[key]
	
	

#################################################
# k-merize reads

namebase=os.path.splitext(ReadsFileName)[0]

statoutPath=('{}_{}_statistics.txt'.format(namebase, ksize))
outPath=('{}_{}_out.txt'.format(namebase, ksize))

statout=open(statoutPath,'w')
out=open(outPath,'w')

statout.write('INPUT\n\nReads file: {}\nReference file: {}\nk-mer size: {}.\n\n'.format(ReadsFileName, RefFileName, ksize))

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


if args.m:
	statout.write("Masking low complexity regions.\n\n")
	remove_lowsh(Hdict)
else:
	statout.write("Not using low complexity mask.\n\n")


plot(Hdict.values())

kmerInfo(Hdict.values())

#print(hash)

###############################################################
# search reference with k-mers, each match gets k-mer value added to the sequence


RefFile=open(RefFileName,'r')

refseqs={}

# printing header for output table:
out.write("Sequence description\tNumber of mapped k-mers\tApproximate number of mapped reads")

nkstat={}
nmapstat={}

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

	nk=refseqs[header]
	nmap=refseqs[header]/ksize

	nkstat[header]=nk
	nmapstat[header]=nmap

	out.write(header + "\t" + str(nk) + "\t" + str(nmap) )

	
nreferences=len(nkstat.keys())
kmean=np.mean(nkstat.values())
mapmean=np.mean(nmapstat.values())

statout.write("\nNumber of reference sequences: {}\nMean number of mapped k-mers per sequence: {}\nMean number of mapped reads per sequence: {}\n".format(nreferences, 
kmean, mapmean))

RefFile.close()
statout.close()
out.close

# TODO:
# has to accept fastq files for reads, not fasta
# test different kmer sizes, how they affect mapping



