#!/usr/bin/env python

import sys

'''
# step 1: open reads file:

arguments = sys.argv
print(arguments)

ReadsFileName=sys.argv[1]

ReadsFile=open(ReadsFileName,'r')


# step 2: k-merize reads

hash={}

for Line in ReadsFile:
	Line=Line.strip() #remove whitespace in beggining and end, because print again adds line endings
	if not '>' in Line:
		length=len(Line)
		length=length-15
		print(length)
		for num in range(0,length):
			first=num
			last=num+15
			kmer=Line[first:last]
			print(kmer)
			hash[kmer]=hash.get(kmer,0)+1
			

print(hash)

ReadsFile.close()
'''
# step 3: search reference with k-mers, each match gets k-mer value added to the sequence
# open reference
# remove newlines in fasta
RefFileName=sys.argv[2]

RefFile=open(RefFileName,'r')

for Line in RefFile:
	Line=Line.strip() #remove whitespace in beggining and end, because print again adds line endings
	print(Line)
RefFile.close()







