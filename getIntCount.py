import argparse
import sys
import HTSeq
import re

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Counts the number of interactions above a below a given distance')
parser.add_argument('-a',type=str,dest="aFile",help="BED file. Output of Stitch-seq_Aligner, where the left and right ends correspond to DNA and RNA, respectively")
parser.add_argument('-d',type=int,dest="distance",help="INT. Minimum distance (in nucleotides) between RNA-DNA ends to be considered as a long range interaccion. Default = 4000 b",default=4000)

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
aFile = args.aFile
distance  = args.distance


##############################################################
# Functions

def getGenomicIntervals(line,distance):
	line=line.strip("\n")
	fields=line.split("\t")
	iv_dna= HTSeq.GenomicInterval(fields[0],int(fields[1]),int(fields[2]),".")
	iv_rna= HTSeq.GenomicInterval(fields[9],int(fields[10]),int(fields[11]),".")
	longRange = False
	if distance ==0:
		longRange = False
	elif iv_dna.chrom != iv_rna.chrom:
		longRange = True
	elif (iv_dna.start - iv_rna.end)>=distance or (iv_rna.start - iv_dna.end)>=distance:
		longRange = True

	return longRange

def getCounts(aFile,distance):
	# Count close and long range interactions

	countCR,countLR = 0,0
	for line in open(aFile,"r"):
		longRange = getGenomicIntervals(line,distance)
		
		if longRange==True:
			countLR += 1
		else:
			countCR += 1

	return [countCR,countLR]

####################################################
# Excecution

print "Counting close and long range interactions per gene"
countCR,countLR = getCounts(aFile,distance) 
print str(countCR)+"\t"+str(countLR)
