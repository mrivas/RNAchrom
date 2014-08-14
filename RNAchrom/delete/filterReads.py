import argparse
import sys
import pysam
import HTSeq
import pysam
import numpy
import re

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Filter reads that are closer than a given distance (interchromosomal pairs are always considered far apart)')
parser.add_argument('-a',type=str,dest="aFile",help="BED file. Output of Stitch-seq_Aligner, where the left and right ends correspond to DNA and RNA, respectively")
parser.add_argument('-b',type=str,dest="bFile",help="BAM file. BAM file to be filter.")
parser.add_argument('-d',type=int,dest="distance",help="INT. Minimum distance (in nucleotides) between RNA-DNA ends to be considered as a long range interaccion. Default = 4000 b",default=4000)
parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file. Default = \"output\"",default="output")

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
aFile = args.aFile
bFile = args.bFile
distance  = args.distance
oFile = args.oFile


##############################################################
# Functions

def getGenomicIntervals(line,distance):
	line=line.strip("\n")
	fields=line.split("\t")
	iv_dna= HTSeq.GenomicInterval(fields[0],int(fields[1]),int(fields[2]),".")
	iv_rna= HTSeq.GenomicInterval(fields[9],int(fields[10]),int(fields[11]),".")
	readName=fields[8]

	longRange = False
	
	# if distance ==0, consider all reads for the analysis=> longRange as True
	if distance == 0:
		longRange = True
	elif iv_dna.chrom != iv_rna.chrom:
		longRange = True
	elif (iv_dna.start - iv_rna.end)>=distance or (iv_rna.start - iv_dna.end)>=distance:
		longRange = True

	return [iv_dna,iv_rna,longRange,readName]

def getKeepReads(aFile,distance):
	keepReads={}
	for line in open(aFile,"r"):
		iv_dna, iv_rna, longRange, readName = getGenomicIntervals(line,distance)
		
		if longRange==False: continue

		keepReads[readName]=1
	
	return keepReads

####################################################
# Excecution


keepReads=getKeepReads(aFile,distance)

bamFile=pysam.Samfile(bFile,"rb")
out=pysam.Samfile(oFile,"wh",header=bamFile.header)

for almnt in bamFile:
	if almnt.qname in keepReads:
		out.write(almnt)

out.close()
