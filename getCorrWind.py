import argparse
import sys
import HTSeq
import numpy
import random
import scipy.stats
import os.path
random.seed(1)

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Computes correlation between number of aimer-target links and number of ChIP-seq/DNase peaks')
parser.add_argument('-l',type=str,dest="lFile",help="BED file. Paired-end alignment reads are presente on each line ")
parser.add_argument('-c',type=str,dest="cFile",help="TAB file. Chromosome size information")
parser.add_argument('-f',type=str,dest="folder",help="STR. Folder where the ChIP-seq/DNase peak files are stored")
parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file. Default = \"output\"",default="output")

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
linksFile = args.lFile
chromFile = args.cFile
folder = args.folder
oFile = args.oFile

windSize=10000
RNAonExon=True
#linksFile="/data2/yu68/bharat-interaction/TwoStep1219-2014-3-17/split_partner/TwoStep1219_fragment_paired_align.txt"
#folder="/data2/rivasas2/rnadna/contrast/ChIP-seq/"
#chromFile="/data2/rivasas2/rnadna/TwoSteps/proofOfConcept/withOverl/exons_only/rmdup/windows/mm9.chrom.sizes"

#################################
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

def linkToIv(line,chromLength,windSize,Type):
# Create iv for DNA and RNA mates 
	line=line.split("\t")

	if Type=="links":
		dna_iv = ivReg( HTSeq.GenomicInterval( line[0],int(line[1]),int(line[2]) ), chromLength, windSize )
		rna_iv = ivReg( HTSeq.GenomicInterval( line[9],int(line[10]),int(line[11]) ), chromLength, windSize )
		if line[16]=="exon" or line[16]=="utr3" or line[16]=="utr5": 
			coding_on=True
		else:
			coding_on=False
		return [dna_iv,rna_iv,coding_on]
	
	else:
		nucl_iv = ivReg( HTSeq.GenomicInterval( line[0],int(line[1]),int(line[2]) ), chromLength, windSize )
		return nucl_iv

def bed2GenomicArrayOfSets(File,chromLength,windSize,RNAonExon,Type):
# Convert DNA-RNA links from BED to GenomicArrayOfSets format

	coverage = {}
	
	for line in open(File,"r"):
		line = line.strip()	
		if Type == "links":
			dna_iv, rna_iv, coding_on = linkToIv(line,chromLength,windSize,Type)
			# Cosiders only RNA mates sitting on exons?
			if RNAonExon and (not coding_on): continue
			if rna_iv in coverage:
				if dna_iv in coverage[rna_iv]:
					coverage[rna_iv][dna_iv] += 1
				else:
					coverage[rna_iv][dna_iv] = 1
			else:
				coverage[rna_iv] = {}
				coverage[rna_iv][dna_iv] = 1
		else:
			nucl_iv = linkToIv(line,chromLength,windSize,Type)
			if nucl_iv in coverage:
				coverage[nucl_iv] += 1
			else:
				coverage[nucl_iv] = 1

	return coverage

##################################
# Execution

chromLength = getchromLength(chromFile)
links = bed2GenomicArrayOfSets(linksFile,chromLength,windSize,RNAonExon,"links") # rna_iv: set of dna_iv
peaks=[]
for fileName in os.listdir(folder):
	if fileName=="procedure.sh": continue
	peaksFile=folder+fileName
	peaks.append( bed2GenomicArrayOfSets(peaksFile,chromLength,windSize,True,"peaks") )# nucleosome_iv: set of peaks_iv

counts={}
countsMatrix={}
# Count hits on targets of region_iv
for chrom in chromLength:
	length=chromLength[chrom]
	for start in range(0,length,windSize):
		end = min(start+windSize,length)

		rna_iv = HTSeq.GenomicInterval(chrom,start,end)
		
		if not rna_iv in links: continue
		if len(links[rna_iv])<5: continue	
		counts[rna_iv]={}
		for dna_iv in links[rna_iv]:
			numPeaks=[]
			for peak in peaks:
				if dna_iv in peak:
					numPeaks.append( peak[dna_iv] )
				else:
					numPeaks.append( 0 )
			
			counts[rna_iv][dna_iv]=[ links[rna_iv][dna_iv], numPeaks]
			if rna_iv in countsMatrix:	
				countsMatrix[rna_iv]=numpy.vstack([countsMatrix[rna_iv],[links[rna_iv][dna_iv]]+ numPeaks])
			else:
				countsMatrix[rna_iv]=[links[rna_iv][dna_iv]]+ numPeaks

out=open(oFile,'w')
for iv in countsMatrix:
	ivString=iv.chrom+":"+str(iv.start)+"-"+str(iv.end)
	output = [ivString] + [ scipy.stats.spearmanr(countsMatrix[iv][:,0],countsMatrix[iv][:,i])[0] for i in range(1,countsMatrix[iv].shape[1]) ]			
	print >>out, "\t".join(map(str,output))
out.close()

# Count histone peaks on region_iv
# Only counts peaks on regions that already have dna counts
# This since most genomic regions having histone peaks aren't
# expected to have interaction with all aimers.				
