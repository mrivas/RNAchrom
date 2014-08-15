#!/usr/bin/env python
import argparse, sys, HTSeq

###########################################################
# Get command line arguments
def getParser():
	parser = argparse.ArgumentParser(description='Get interesting read\'s IDs (printed to STDOUT). That\'s, reads that: (1) have spliced junctions concordant with the annotated genome, and (2) have a mate farther apart than a given distance')
	parser.add_argument('-b',type=str,dest="bFile",help="BAM file(s). A single aligment file contains all paired-end mates.")
	parser.add_argument('-d',type=int,dest="dist",help="INT. Minimum distance between paired-end mates")
	parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file.")

	if len(sys.argv) == 1:
		print >> sys.stderr,parser.print_help()
		exit(0)
	return parser
#############################################################
def main():
	args = getParser().parse_args()
	bFile = args.bFile
	dist = args.dist
	oFile = args.oFile

	out=open(oFile,'w')
	first = None

	bam=HTSeq.BAM_Reader( bFile )
	# Filter by regions of interest
	for almnt in bam:

		# Discard not aligned reads, multihits
		if (not almnt.aligned) or almnt.optional_field("NH")>1: continue
		if first == None:    
			first = almnt
		else: 
			second = almnt
			# Check that first and second reads have same read name
			if first.read.name != second.read.name:
				first=second
				continue
			
			# Skip self-ligating mates
			if abs(first.inferred_insert_size) <= dist: 
				first=None
				continue
			
			sjFirst  = first.optional_field("jM")
			sjSecond = second.optional_field("jM")
			
			printRead=False
			for sj in sjFirst:
				if sj>=20: printRead=True
			for sj in sjSecond:
				if sj>=20: printRead=True

			# Print mates if at least one of them has a database spliced-junction 
			if printRead: 
				print >>out,first.read.name

			first = None
	out.close()
###############################################################
if __name__ == '__main__':
	main()
