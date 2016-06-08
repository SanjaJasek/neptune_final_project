#!/usr/bin/env python

import sys

# step 1: open file:

arguments = sys.argv
print(arguments)

InFileName=sys.argv[1]

InFile=open(InFileName,'r')


# step 2: k-merize reads

hash={}

for Line in InFile:
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


InFile.close()
