#!/usr/bin/env python
import HTSeq, argparse, sys
import RNAchrom

def getParser():
	# Get command line arguments
	parser = argparse.ArgumentParser(description='Find strong aimer-target relations. It uses aware-links (know and inferred) as seed to determine target genomic regions(TGR). A TGR is build by overlapping the windows of target DNA-mates. Conversely, exons are used as aimer genomic regions')
	parser.add_argument('-a',type=str,dest="aFile",help="TAB file. TAB separated file with aware mates' information (output of annotateBAM.py).")
	parser.add_argument('-b',type=str,dest="bFile",help="TAB file. TAB separated file with blind mates' information (output of annotateBAM.py).")
	parser.add_argument('-g',type=str,dest="gFile",help="GTF file.")
	parser.add_argument('-w',type=int,dest="windSize",help="INT. Window size (nt) to build around aware links to determine target regions.")
	parser.add_argument('-d',type=int,dest="dist",help="INT. Distance between mates used as lower threshold to call them non self-ligating.")
	parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file")

	if len(sys.argv) == 1:
		print >> sys.stderr,parser.print_help()
		exit(0)
	return parser

###########################################################
def main():
#	args = getParser().parse_args()
#	aFile = args.aFile
#	bFile = args.bFile
#	gtfFile = args.gFile
#	windSize = args.windSize
#	dist = args.dist
#	oFile = args.oFile

	aFile="/data2/rivasas2/rnadna/TwoSteps/allReads/annotations/mm9/TwoStep1024.aware_rmdup_annotatedBAM.bed"
	bFile="/data2/rivasas2/rnadna/TwoSteps/allReads/annotations/mm9/TwoStep1024.blind_rmdup_annotatedBAM.bed"
	gtfFile="/home/rivasas2/tools/genomes/mouse/mm9/Mus_musculus.NCBIM37.67_chr.gtf"
	windSize=10000
	dist=10000
	oFile="test.txt"

	# Classify mates as aware or blind. Ambiguous mates are discarded
	# Aware mates are: know and infered. However, if the DNA mates has a know sj the aware mates is dicarded
	# Known aware mates have a RG saying so
	# Inferred mates are the ones when only one mate have a known sj (jM>20 on STAR)
	# Ambiguous mates are those where both mates show SJs or non of them overlap exons.
	# Self-ligating mates are also discarded
	print "Loading known aware links"
	awareMates = RNAchrom.loadMates(aFile,"aware",dist)
	print "Loading blind links"
	inferredMates,blindMates = RNAchrom.loadMates(bFile,"blind",dist)
	print "Merging known and inferred aware links"
	awareMates.update(inferredMates)
	print "Creating genes' IV"
	genesID_IV,geneIV_ID=RNAchrom.genes(gtfFile)
	print "Counting and printing results"
	RNAchrom.buildLinks(awareMates,blindMates,windSize,genesID_IV,oFile)	
####################################################
if __name__ == '__main__':
	main()
