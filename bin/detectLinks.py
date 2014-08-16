#!/usr/bin/env python
import HTSeq, argparse, sys
import RNAchrom

def getParser():
	# Get command line arguments
	parser = argparse.ArgumentParser(description='Find strong aimer-target relations. It uses aware-links (know and inferred) as seed to determine target genomic regions(TGR). A TGR is build by overlapping the windows of target DNA-mates. Conversely, exons are used as aimer genomic regions')
parser.add_argument('-m',type=str,dest="mFile",help="TAB file. TAB separated file with mates' information (output of annotateBAM2.py).")
	parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file. Default = \"output\"",default="output")

	if len(sys.argv) == 1:
		print >> sys.stderr,parser.print_help()
		exit(0)
	return parser

###########################################################
def main():
	args = getParser().parse_args()
	mFile = args.mFile
	oFile = args.oFile
	# Classify mates as aware or blind. Ambiguous mates are discarded
	# Aware mates are: know and infered. However, if the DNA mates has a know sj the aware mates is dicarded
	# Known aware mates have a RG saying so
	# Inferred mates are the ones when only one mate have a known sj (jM>20 on STAR)
	# Ambiguous mates are those where both mates show SJs or non of them overlap exons.
	# Self-ligating mates are also discarded
	mates = RNAchrom.loadMates(mFile)
	# Build links
	links = RNAchrom.buildLinks(mates)	

	# Print results
	out = open(oFile, 'w')
	for geneID in links:
		for link in links[geneID]:
			print >>out, link.printFormat()
	out.close()

mates[ readName ] # mate indexed by read name 
	Type # either aware, blind, or ambiguous
	mate1.iv # 
	mate1.annot # itron/exon/.|repeat/.|sj/. 
	mate2.iv # 
	mate2.annot # itron/exon/.|repeat/.|sj/. 

links[geneID] = [link1, link2, ...]
for link in links[geneID]: #fields 
	link['aimerIv']     # aimer iv including intronic regions
	link['aimerLength'] # length coding regions only
	link['aimerAnnot']   # genes overlapping aimer.iv
	link['targetIv'] # iv target region
	link['targetLength'] # length of target region
	link['targetAnnot'] # repeats and genes overlapping target.iv
	link['targetAwareCounts'] # number of aware-links hitting target
	link['targetBlindCounts'] # number of blind-links hitting target 
	link['targetSEM'] # standard error of the mean (positions of the target-hits)



####################################################
if __name__ == '__main__':
	main()
