import argparse
import sys
import HTSeq
import numpy
import random
import scipy.stats
import os.path
import pickle
random.seed(1)

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Computes correlation between number of aimer-target links and number of ChIP-seq/DNase peaks')
parser.add_argument('-l',type=str,dest="lFile",help="BED file. Paired-end alignment reads are presente on each line. Output of annotateBAM.py ")
parser.add_argument('-e',type=str,dest="exon",help="STR BOOL. \"True\", for checking that a mate should overlap and exon to be counted as outgoing link of a window (if both mates are on the same window, at least one of them should overlap an exon). \"False\", if want to count all links.")
parser.add_argument('-w',type=int,dest="window",help="INT. Window size (nt) when splitting genome.")
parser.add_argument('-d',type=int,dest="dist",help="INT. Distance between mates to be consider not selfligating.")
parser.add_argument('-c',type=str,dest="cFile",help="TAB file. Chromosome size information")
parser.add_argument('-f',type=str,dest="folder",help="STR. Folder where the ChIP-seq/DNase peak files are stored")
parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file." )

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
linksFile = args.lFile
checkExon = args.exon
windSize = args.window
dist = args.dist
chromFile = args.cFile
folder = args.folder
oFile = args.oFile

#windSize=10000
#checkExon=False
#dist=2000
#linksFile="/data2/rivasas2/rnadna/TwoSteps/allReads/TwoStep1024.annotatedBAM.bed"
#folder="/data2/rivasas2/rnadna/contrast/ChIP-seq/"
#chromFile="/data2/rivasas2/rnadna/TwoSteps/proofOfConcept/withOverl/exons_only/rmdup/windows/mm9.chrom.sizes"
#oFile="test.bed"

##################################
# Functions

def getchromLength(chromFile):
	chromLength={}
	for line in open(chromFile,'r'):
		line=line.strip().split("\t")
		chromLength[line[0]] = int(line[1])
	return chromLength

def ivReg(iv,chromLength,windSize):
# Regularization of intervals. The iv are transformed 
# to canonical regions. This to avoid partial overlappings, and to get 
# well defined aimer and target windows
	chrom = iv.chrom
	start = iv.start
	end   = iv.end
	if chrom in chromLength:
		length = chromLength[chrom]
	
		start_s = start - start % windSize
		start_e = min(start_s + windSize, length)
		end_s   = end - end % windSize
		end_e   = min(end_s   + windSize, length)
		
		# Return the canonical interval with the largest overlap,
		# or chose one at random if iv is equally present on two of them
		if (start_e-start) > (end-end_s):
			reg_iv = HTSeq.GenomicInterval(chrom,start_s,start_e)	
		elif (start_e-start) < (end-end_s):
			reg_iv = HTSeq.GenomicInterval(chrom,end_s,end_e)	
		else:
			if random.random()>0.5:
				reg_iv = HTSeq.GenomicInterval(chrom,start_s,start_e)	
			else:
				reg_iv = HTSeq.GenomicInterval(chrom,end_s,end_e)	
		
		return reg_iv
	else: # Current read not on known chromosome
		return "NoPresent"

def bed2Peaks(File,chromLength,windSize):
# Convert DNA-RNA links from BED to GenomicArrayOfSets format

	coverage = {}
	
	for line in open(File,"r"):
		line = line.strip().split("\t")	
		iv = ivReg( HTSeq.GenomicInterval( line[0],int(line[1]),int(line[2]) ), chromLength, windSize )
		if type(iv)!=HTSeq._HTSeq.GenomicInterval: continue
		if iv in coverage:
			coverage[iv] += 1
		else:
			coverage[iv] = 1

	return coverage

def lineToIv(line,chromLength,windSize,dist):
# Create iv for DNA and RNA mates 
	line=line.split("\t")
	coding1,coding2=False,False
	selfLig=True
	iv1 = ivReg( HTSeq.GenomicInterval( line[0],int(line[1]),int(line[2]) ), chromLength, windSize )
	iv2 = ivReg( HTSeq.GenomicInterval( line[5],int(line[6]),int(line[7]) ), chromLength, windSize )
	if type(iv1)!=HTSeq._HTSeq.GenomicInterval or type(iv2)!=HTSeq._HTSeq.GenomicInterval:
		return ["NoPresent","NoPresent","NoPresent","NoPresent","NoPresent"]
	else:
		# check if mates overlap exons
		if line[4].split("|")[-2]=="exon": 
			coding1=True
		if line[9].split("|")[-2]=="exon": 
			coding2=True
		# Check if mates are selfligating
		if iv1.chrom!=iv2.chrom:
			selfLig=False
		elif abs(iv1.start-iv2.start)>dist:
			selfLig=False

		return [iv1,iv2,coding1,coding2,selfLig]
	

def bed2Links(bFile,chromLength,windSize,checkExon,dist):
# Convert DNA-RNA links from BED to GenomicArrayOfSets format

	links = {}
	
	for line in open(bFile,"r"):
		line = line.strip()	
		iv1, iv2, coding1,coding2,selfLig = lineToIv(line,chromLength,windSize,dist)
		# Ignore links not present on known chromosomes
		if type(iv1)!=HTSeq._HTSeq.GenomicInterval or type(iv2)!=HTSeq._HTSeq.GenomicInterval: continue
		# Ignore selfLigating links (mates closer to each other by less than 2k nt)
		if selfLig: continue

		if iv1!=iv2:
			# Count only if iv1 overlap an exon
			if (not checkExon) or coding1:
				if iv1 in links:
					if iv2 in links[iv1]:
						links[iv1][iv2] += 1
					else:
						links[iv1][iv2] = 1
				else:
					links[iv1] = {}
					links[iv1][iv2] = 1
			# Count only if iv2 overlap an exon
			if (not checkExon) or coding2:
				if iv2 in links:
					if iv1 in links[iv2]:
						links[iv2][iv1] += 1
					else:
						links[iv2][iv1] = 1
				else:
					links[iv2] = {}
					links[iv2][iv1] = 1
		else:
			# If iv1==iv2 at least one end must overlap and exon
			if (not checkExon) or ( coding1 or coding2 ):
				if iv1 in links:
					if iv2 in links[iv1]:
						links[iv1][iv2] += 1
					else:
						links[iv1][iv2] = 1
				else:
					links[iv1] = {}
					links[iv1][iv2] = 1
			

	return links

##################################
# Execution

print "Loading chromosomes"
chromLength = getchromLength(chromFile)
print "Creating links database"
links = bed2Links(linksFile,chromLength,windSize,checkExon,dist) # rna_iv: set of dna_iv
print "Creating peaks database"
peaks=[]
header=["coord"]
for fileName in os.listdir(folder):
	if fileName=="procedure.sh": continue
	header.append(fileName)
	peaksFile=folder+fileName
	peaks.append( bed2Peaks(peaksFile,chromLength,windSize) )# nucleosome_iv: set of peaks_iv

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
