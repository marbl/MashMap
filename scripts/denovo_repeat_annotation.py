################################################################################
# @file     denovo_repeat_annotation.py
# @author   Chirag Jain (chirag.jain@nih.gov)
# @purpose  postprocess mashmap output for denovo repeat annotation
#  
#						A common task when analysing repeats in a genome is to label
#						repetitive segments in the genome, i.e., segments that are similar
#						to other loci. Suppose you want to annotate repeats of 
#						mimimum length X and minimum identity Y%, mashmap (and this script)
#						can provide a bed file specifying repeats in the genome

#						Step 1: Run mashmap for self-alignment of genome 
#						Step 2: Run this script to obtain a bed file
#						(optional) Step 3: Merge duplicate intervals using bedtools

#						Command for each step (suppose you wish X=5000, Y=95%):

#						mashmap -r genome.fa -q genome.fa -f none -s 5000 --pi 95
#						python denovo_repeat_annotation.py mashmap.out 5000 95 > repeats.bed
#						bedtools merge -i repeats.bed > repeats.merged.bed

#						Note: This is a fast approximation using an alignment-free method
################################################################################

from os import sys

CHROMOSOMECOL1=0
STARTCOL1=2
ENDCOL1=3
STRAND=4
CHROMOSOMECOL2=5
STARTCOL2=7
ENDCOL2=8
IDENTITY=9

repeatList = []

with open(sys.argv[1]) as f:
	for line in f:
		rowElements = line.split()
		chromosome1 = rowElements[CHROMOSOMECOL1]
		start1 = int(rowElements[STARTCOL1])
		end1 = int(rowElements[ENDCOL1])
		strand = rowElements[STRAND]
		chromosome2 = rowElements[CHROMOSOMECOL2]
		start2 = int(rowElements[STARTCOL2])
		end2 = int(rowElements[ENDCOL2])
		identity = float(rowElements[IDENTITY])

		if (chromosome1 != chromosome2 or (abs(start1 - start2) >= 1.5 * int(sys.argv[2]) and abs(end1 - end2) >= 1.5 * int(sys.argv[2]))):
			if (end1 - start1 + 1 >= int(sys.argv[2]) and identity + 1 >= float(sys.argv[3])):	#added one to identity for sensitivity
				repeatList.append((chromosome1, start1, end1))

# sort the intervals
repeatList.sort()

# print the intervals
for (chromosome, start, end) in repeatList:
	print(chromosome +  "\t" + str(start) + "\t" +  str(end+1))

# Users can optionally merge the output using bedtools (see instructions at the top of this file)
