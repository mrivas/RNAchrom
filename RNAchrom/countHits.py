#!/usr/bin/env python
import argparse, sys, HTSeq
import RNAchrom

###########################################################
# Get command line arguments
def getParser():
	parser = argparse.ArgumentParser(description='Counts the number of RNA and DNA mates per gene (including introns). Only mates where the RNA-ends overlap exons are use.')
	parser.add_argument('-a',type=str,dest="aFile",help="BED file. Output of annotateBAM.py, where the left and right ends correspond to DNA and RNA, respectively")
	parser.add_argument('-g',type=str,dest="gFile",help="GTF file. Annotation file with biotype information")
	parser.add_argument('-d',type=int,dest="distance",help="INT. Minimum distance (in nucleotides) between RNA-DNA ends to be considered as a long range interaccion. Default = 2000 b",default=2000)
	parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file.")

	if len(sys.argv) == 1:
		print >> sys.stderr,parser.print_help()
		exit(0)
	return parser
#########################################################
def main():
	args = getParser().parse_args()
	aFile = args.aFile
	gtfFile = args.gFile
	distance  = args.distance
	oFile = args.oFile

	#aFile="/data2/yu68/bharat-interaction/TwoStep1024-2014-3-14/split_partner/TwoStep1024_fragment_paired_align.txt"
	#gtfFile="/home/rivasas2/tools/genomes/mouse/mm9/Mus_musculus.NCBIM37.67.gtf.gz"
	#distance=5000
	#oFile="test.txt"
	print "Getting genes"
	geneID_IV,geneIV_ID = RNAchrom.genes(gtfFile)
	print "Counting long range interactions per gene"
	countsDNA, countsRNA = RNAchrom.getCounts(aFile,distance,geneIV_ID) 
	print "Getting biotypes"
	bioType = RNAchrom.getBioType(gtfFile)
	print "Saving results to: "+oFile
	out = open(oFile,'w')
	for geneID in geneID_IV:
		length = geneID_IV[geneID].end-geneID_IV[geneID].start
		coord = geneID_IV[geneID].chrom+":"+str(geneID_IV[geneID].start)+"-"+str(geneID_IV[geneID].end)+";"+geneID_IV[geneID].strand

		if geneID in countsDNA: 
			cDNA=str(countsDNA[geneID])
		else:
			cDNA = "0"
		if geneID in countsRNA: 
			cRNA=str(countsRNA[geneID])
		else:
			cRNA = "0"
		
		print >>out,  geneID+"\t"+cDNA+"\t"+cRNA+"\t"+str(length)+"\t"+bioType[geneID]+"\t"+coord

	out.close()
##########################################################
if __name__ == '__main__':
	main()
