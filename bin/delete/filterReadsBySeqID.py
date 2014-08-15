import argparse
import sys
import pysam
import HTSeq

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Filter out sequence IDs.')
parser.add_argument('-i',type=str,dest="iFile",help="BED file. Coordenate-sorted output of Stitch-seq_Aligner.")
parser.add_argument('-s',type=str,dest="sFile",help="TXT file. List of sequence IDs to be keept.")
parser.add_argument('-o',type=str,dest="oFile",help="BED file. Name of output file. Same format as Stitch-seq_Aligner output")

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
iFile = args.iFile
sFile = args.sFile
oFile = args.oFile

#####################################################################################3
# Execution

##############
# Create hash of sequence IDs.
seqID={}
for line in open(sFile,'r'):
	line=line.strip()
	seqID[line]=1

#############
# Check if seq ID is on Stitch-seq file
out=open(oFile,'w')

for line in open(iFile,"r"):
	line=line.strip()
	ID=line.split("\t")[8]
	if ID in seqID:
		print >>out, line

out.close()
