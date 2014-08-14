import argparse
import sys
import pysam
import HTSeq

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Remove overlapping reads from Stitch-seq_Aligner output.')
parser.add_argument('-i',type=str,dest="iFile",help="BED file. Coordenate-sorted output of Stitch-seq_Aligner.")
parser.add_argument('-o',type=str,dest="oFile",help="BED file. Name of output file. Same format as Stitch-seq_Aligner output")

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
iFile = args.iFile
oFile = args.oFile

#####################################################################################3
# Execution

out=open(oFile,'w')

firstIVa,firstIVb = None,None

for line in open(iFile,"r"):
	line=line.strip("\n")
	fields=line.split("\t")

	if firstIVa==None:
		firstIVa = HTSeq.GenomicInterval(fields[0],int(fields[1]),int(fields[2]),fields[3] )
		firstIVb = HTSeq.GenomicInterval(fields[9],int(fields[10]),int(fields[11]),fields[12] )
		firstLine = line
	else:
		secondIVa = HTSeq.GenomicInterval(fields[0],int(fields[1]),int(fields[2]),fields[3] )
		secondIVb = HTSeq.GenomicInterval(fields[9],int(fields[10]),int(fields[11]),fields[12] )
		
		printed=False
		if secondIVa!=firstIVa or secondIVb!=firstIVb:
			print >>out, firstLine
			printed=True
		
		firstIVa = secondIVa
		firstIVb = secondIVb
		firstLine = line
	
if printed==True: print >>out, firstLine
out.close()
