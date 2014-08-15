#!/usr/bin/env python
import argparse, HTSeq, numpy, scipy.stats, pickle, random, sys, os
import RNAchrom
random.seed(1)

###########################################################
# Get command line arguments
def getParser():
	parser = argparse.ArgumentParser(description='Computes correlation between number of aimer-target links and number of ChIP-seq/DNase peaks')
	parser.add_argument('-l',type=str,dest="lFile",help="BED file. Paired-end alignment reads are presente on each line. Output of annotateBAM.py ")
	parser.add_argument('-e',type=str,dest="exon",help="STR BOOL. \"True\", for checking that a mate should overlap and exon to be counted as outgoing link of a window (if both mates are on the same window, at least one of them should overlap an exon). \"False\", if want to count all links.")
	parser.add_argument('-w',type=int,dest="window",help="INT. Window size (nt) when splitting genome.")
	parser.add_argument('-t',type=str,dest="linkType",help="STR. Link type: aware, blind.")
	parser.add_argument('-d',type=int,dest="dist",help="INT. Distance between mates to be consider not selfligating.")
	parser.add_argument('-c',type=str,dest="cFile",help="TAB file. Chromosome size information")
	parser.add_argument('-f',type=str,dest="folder",help="STR. Folder where the ChIP-seq/DNase peak files are stored")
	parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file." )

	if len(sys.argv) == 1:
		print >> sys.stderr,parser.print_help()
		exit(0)
	return parser

###########################################################
def main():
	# Arguments
	args = getParser().parse_args()
	linksFile = args.lFile
	checkExon = args.exon
	if checkExon=="True": checkExon=True
	else:                 checkExon=False
	windSize = args.window
	dist = args.dist
	linkType = args.linkType
	chromFile = args.cFile
	folder = args.folder
	oFile = args.oFile

	print "Loading chromosomes"
	chromLength = RNAchrom.chromLength(chromFile)
	print "Creating links database"
	links = RNAchrom.bed2links(linksFile,chromLength,windSize,checkExon,dist,linkType) # rna_iv: set of dna_iv
	print "Creating peaks database"
	peaks=[]
	header=["coord"]
	for fileName in os.listdir(folder):
		if fileName=="procedure.sh": continue
		header.append(fileName)
		peaksFile=folder+fileName
		peaks.append( bed2peaks(peaksFile,chromLength,windSize) )# nucleosome_iv: set of peaks_iv

	print "Matching number of hits and peaks"
	countsMatrix={}
	# Count hits on targets of region_iv
	# Count histone peaks on region_iv
	# Only counts peaks on regions that already have links counts
	# This since most genomic regions having histone peaks aren't
	# expected to have interaction with all aimers.				
	for iv1 in links: # Iterate over aimer regions
		# Requieres that current aimer has at least 5 targets
		if len(links[iv1])<5: continue
		for iv2 in links[iv1]: # Iterate over targets of current aimer
			numPeaks=[]
			for peak in peaks: # Fetch histone information
				# Fetch of histone peaks for current target region
				if iv2 in peak:
					numPeaks.append( peak[iv2] )
				else:
					numPeaks.append(0)
			# For current amimer, saves number of hits (col 1) and number of peaks for each target (rows)	
			if iv1 in countsMatrix:
				countsMatrix[iv1]=numpy.vstack( [ countsMatrix[iv1], [links[iv1][iv2]]+ numPeaks] )
			else:
				countsMatrix[iv1]=[ links[iv1][iv2]  ] + numPeaks

	print "Printing results to output file"
	out=open(oFile,'w')
	print >>out, "\t".join(map(str,header))
	for iv in countsMatrix:
		ivString=iv.chrom+":"+str(iv.start)+"-"+str(iv.end)
		output = [ivString] + [ scipy.stats.spearmanr(countsMatrix[iv][:,0],countsMatrix[iv][:,i]) for i in range(1,countsMatrix[iv].shape[1]) ]			
		print >>out, "\t".join(map(str,output))
	out.close()

	print "Save objects to file"
	out=open(oFile+"_countsMatrix.pkl","w")
	pickle.dump(countsMatrix,out)
	out.close()
	out=open(oFile+"_links.pkl","w")
	pickle.dump(links,out)
	out.close()
	out=open(oFile+"_peaks.pkl","w")
	pickle.dump(peaks,out)
	out.close()

##################################################################
if __name__ == "__main__":
	main()

