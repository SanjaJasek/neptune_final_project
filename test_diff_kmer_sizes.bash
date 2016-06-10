#!/bin/bash -x

mapper="./k-mer_mapper.py"
reads="cow_mRNA_reads.fasta"
ref="cow_mRNA.fasta"

# collect statistics for no mask
nomaskHs=()
nomaskks=()
nomaskrs=()
nomaskout="nomask_stats.txt"

# collect statistics for mask
maskHs=()
maskks=()
maskrs=()
maskout="mask_stats.txt"

# write headers:
header () {
fil=$1
echo -e "k-mer size\tmean Shannon entorpy\tmean number of k-mers per sequence\tmean number of mapped reads per sequence" > $fil
}

header $nomaskout
header $maskout

stats () {
stfil=$1
string=$2

value=$( grep "$dtring" "$stfil" )
value=${value##*: }

#return $value
#echo $value
}

for thr in 5 10 15 20 25 30
do 
	# map reads without mask
	nomask="cow_mRNA_reads_${thr}_statistics.txt"

	[[ -e $nomask ]] || $mapper -reads $reads -ref $ref -k $thr

	# get stats
	
	#nomaskH=$(stats $nomask "Shannon entropy")
	#nomaskHs+=($nomaskH)
	stats $nomask "Shannon entropy"
	nomaskH=$value
	
	#nomaskk=$( stats $nomask "k-mers per sequence"  )
	#nomaskks+=($nomaskk)
	stats $nomask "k-mers per sequence" 
	nomaskk=$value

	#nomaskr=$( stats $nomask "reads per sequence"  )
	#nomaskrs+=($nomasrk)
	stats $nomask "reads per sequence"
	nomaskr=$value

	echo -e "$nomaskH\t$nomaskk\t$nomaskr" >> $nomaskout

	# map reads with mask

	mask="cow_mRNA_reads_${thr}_mask_statistics.txt"

	[[ -e $mask ]] || $mapper -reads $reads -ref $ref -k $thr -m

	# get stats
	
	#maskH=$(stats $mask "Shannon entropy")
	#maskHs+=($maskH)
	stats $mask "Shannon entropy"
	maskH=$value
	
	#maskk=$( stats $mask "k-mers per sequence"  )
	#maskks+=($maskk)
	stats $mask "k-mers per sequence" 
	maskk=$value

	#maskr=$( stats $mask "reads per sequence"  )
	#masrks+=($masrk)
	stats $mask "reads per sequence"
	maskr=$value

	echo -e "$maskH\t$maskk\t$maskr" >> $maskout
done

