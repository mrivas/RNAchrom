#!/usr/bin/env python
import argparse, sys, HTSeq, pickle, random 
import RNAchrom
random.seed(1)

###########################################################
# Get command line arguments
def getParser():
	parser = argparse.ArgumentParser(description='Computes fold change of blind-links on regions targeted by aware-lings vs blind-links targeting everywhere')
	parser.add_argument('-a',type=str,dest="aFile",help="BED file. Aware-links file. The first and second mates correspond to DNA and RNA annotated using annotateBAM.py ")
	parser.add_argument('-b',type=str,dest="bFile",help="BED file. Blind-links file. annotated using annotateBAM.py ")
	parser.add_argument('-e',type=str,dest="exon",help="STR BOOL. \"True\", for checking that a mate should overlap and exon to be counted as outgoing link of a window (if both mates are on the same window, at least one of them should overlap an exon). \"False\", if want to count all links.")
	parser.add_argument('-w',type=int,dest="window",help="INT. Window size (nt) when splitting genome.")
	parser.add_argument('-d',type=int,dest="dist",help="INT. Distance between mates to be consider not selfligating.")
	parser.add_argument('-c',type=str,dest="cFile",help="TAB file. Chromosome size information")
	parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file." )

	if len(sys.argv) == 1:
		print >> sys.stderr,parser.print_help()
		exit(0)
	return parser
##################################
# Execution
def main():

	args = getParser().parse_args()
	awareFile = args.aFile
	blindFile = args.bFile
	checkExon = args.exon
	if checkExon=="True": checkExon=True
	else:                 checkExon=False
	windSize = args.window
	dist = args.dist
	chromFile = args.cFile
	oFile = args.oFile

	if checkExon=="True":
		checkExon=True
	else:
		checkExon=False
	#windSize=10000
	#checkExon="False"
	#dist=2000
	#awareFile="/data2/rivasas2/rnadna/TwoSteps/allReads/correlations/mm9/TwoStep1024.aware_rmdup_annotatedBAM.bed"
	#blindFile="/data2/rivasas2/rnadna/TwoSteps/allReads/correlations/mm9/TwoStep1024.all_rmdup_annotatedBAM.bed"
	#chromFile="/data2/rivasas2/rnadna/TwoSteps/allReads/blindVsAwareLinks/mm9/mm9.chrom.sizes"
	#oFile="test.bed"

	print "Loading chromosomes"
	chromLength = RNAchrom.chromLength(chromFile)
	print "Creating links database"
	awareLinks = RNAchrom.bed2links(awareFile,chromLength,windSize,checkExon,dist,"aware") # rna_iv as iv1
	blindLinks = RNAchrom.bed2links(blindFile,chromLength,windSize,checkExon,dist,"blind") #
	print "Matching number of hits and peaks"
	ratio={}
	awareCounts,specificBlindCounts,allBlindCounts={},{},{}
	# Count hits on targets of region_iv
	for iv1 in awareLinks: # Iterate over aimer regions
		# Count all blind-link targets
		allBlindCount=0
		#if iv1 in blindLinks:
		for iv2 in blindLinks[iv1]:
			allBlindCount += blindLinks[iv1][iv2]

		# Count blind-links on aware-links targets
		specificBlindCount=0
		#if iv1 in blindLinks:
		for iv2 in awareLinks[iv1]: # Iterate over targets of current aimer
			if iv2 in blindLinks[iv1]:
				specificBlindCount += blindLinks[iv1][iv2]
		# Count aware-links targets
		awareCount=0
		for iv2 in awareLinks[iv1]: # Iterate over targets of current aimer
			awareCount += awareLinks[iv1][iv2]
		
		ratio[iv1] = float(specificBlindCount) / float(allBlindCount)		
		
		awareCounts[iv1]        = awareCount
		specificBlindCounts[iv1] = specificBlindCount
		allBlindCounts[iv1]     = allBlindCount

	print "Printing results to output file"
	out=open(oFile,'w')
	for iv in ratio:
		print >>out, iv.chrom+"\t"+str(iv.start)+"\t"+str(iv.end)+"\t"+str(ratio[iv])+"\t"+str(awareCounts[iv])+"\t"+str(specificBlindCounts[iv])+"\t"+str(allBlindCounts[iv])
				
	out.close()

	print "Save objects to file"
	out=open(oFile+"_blindLinks.pkl","w")
	pickle.dump(blindLinks,out)
	out.close()
	out=open(oFile+"_awareLinks.pkl","w")
	pickle.dump(awareLinks,out)
	out.close()
	out=open(oFile+"_ratio.pkl","w")
	pickle.dump(ratio,out)
	out.close()
###############################################################3
if __name__ == '__main__':
	main()

