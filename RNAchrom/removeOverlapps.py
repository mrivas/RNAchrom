#!/usr/bin/env python
import argparse, sys, HTSeq
import RNAchrom

###########################################################
# Get command line arguments
def getParser():
	parser = argparse.ArgumentParser(description='Remove overlapping reads.')
	parser.add_argument('-i',type=str,dest="iFile",help="BED file. Coordenate-sorted output of annotateBAM.py.")
	parser.add_argument('-o',type=str,dest="oFile",help="BED file. Name of output file. Same format as input file.")

	if len(sys.argv) == 1:
		print >> sys.stderr,parser.print_help()
		exit(0)
	return parser
#####################################################################################3
def main():
	args = getParser().parse_args()
	iFile = args.iFile
	oFile = args.oFile

	# Execution

	out=open(oFile,'w')

	firstIVa,firstIVb = None,None

	for line in open(iFile,"r"):
		line=line.strip("\n")
		fields=line.split("\t")

		if firstIVa==None:
			firstIVa = HTSeq.GenomicInterval(fields[0],int(fields[1]),int(fields[2]),fields[3] )
			firstIVb = HTSeq.GenomicInterval(fields[6],int(fields[7]),int(fields[8]),fields[9] )
			firstLine = line
		else:
			secondIVa = HTSeq.GenomicInterval(fields[0],int(fields[1]),int(fields[2]),fields[3] )
			secondIVb = HTSeq.GenomicInterval(fields[6],int(fields[7]),int(fields[8]),fields[9] )
			
			printed=False
			if secondIVa!=firstIVa or secondIVb!=firstIVb:
				print >>out, firstLine
				printed=True
			
			firstIVa = secondIVa
			firstIVb = secondIVb
			firstLine = line
		
	if printed==True: print >>out, firstLine
	out.close()
#####################################################################################3
if __name__ == '__main__':
	main()
