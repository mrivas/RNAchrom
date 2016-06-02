#!/usr/bin/env python
import HTSeq, argparse, sys
import RNAchrom

def getParser():
	# Get command line arguments
	parser = argparse.ArgumentParser(description='Annotate both mates of paired-end alignments')
	parser.add_argument('-b',type=str,dest="bFile",help="BAM file(s). A single aligment file contains all paired-end mates, or a comma separated list of 2 files, one per each set of mates (useful when mate origin is important).")
	parser.add_argument('-g',type=str,dest="gFile",help="GTF file. Annotation file with biotype information")
	parser.add_argument('-r',type=str,dest="rFile",help="BED file. Annotation file with repeats information")
	parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file. Default = \"output\"",default="output")

	if len(sys.argv) == 1:
		print >> sys.stderr,parser.print_help()
		exit(0)
	return parser

###########################################################
def main():
	args = getParser().parse_args()
	bFile = args.bFile
	gFile = args.gFile
	rFile = args.rFile
	oFile = args.oFile
	# Excecution #####################################################
	print "Getting repeats"
	repeat = RNAchrom.repeats(rFile)
	print "Getting genes,biotype, and exons"
	gene,exon = RNAchrom.geneExonIV(gFile)
	print "Saving results to: "+oFile
	out = open(oFile,'w')
	first = None
	firstMate={}
	for i in range( len(bFile.split(",") ) ):
		bam=HTSeq.BAM_Reader( bFile.split(",")[i] )
		
		for almnt in bam:
			###########################################################
			# If one file with both mates
			if len(bFile.split(","))==1:
				# Discard not aligned reads, and multihits
				if (not almnt.aligned) or almnt.optional_field("NH")>1: continue
				if first == None:    
					first = almnt
				else: 
					second = almnt
					# Check that first and second reads have same read name
					if first.read.name != second.read.name:
						first=second
						continue
					annotFirst  = RNAchrom.annotate(first,gene,exon,repeat)
					annotSecond = RNAchrom.annotate(second,gene,exon,repeat)
					
					if first.pe_which=="first": # Always print first in pair at the beginning
						output = RNAchrom.formatOutput(first,second,annotFirst,annotSecond)
					else:
						output = RNAchrom.formatOutput(second,first,annotSecond,annotFirst)
					print >>out,output  
					first = None
			###########################################################
			# If two files, each with one mate
			elif len(bFile.split(","))>1:
				# Discard not aligned reads, and multihits
				if (not almnt.aligned) or almnt.optional_field("NH")>1: continue
				annot = RNAchrom.annotate(almnt,gene,exon,repeat)
				
				if i==0:
					firstMate[almnt.read.name]=[almnt]
					firstMate[almnt.read.name].append( annot )
				elif i>0 and (almnt.read.name in firstMate):
					output = RNAchrom.formatOutput(firstMate[almnt.read.name][0],almnt,firstMate[almnt.read.name][1],annot)
					print >>out, output

	out.close()

####################################################
if __name__ == '__main__':
	main()
