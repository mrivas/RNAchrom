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

parser = argparse.ArgumentParser(description='Computes fold change of blind-links on regions targeted by aware-lings vs blind-links targeting everywhere')
parser.add_argument('-a',type=str,dest="aFile",help="BED file. Aware-links file. The first and second mates correspond to DNA and RNA. ")
parser.add_argument('-b',type=str,dest="bFile",help="BED file. Blind-links file. ")
parser.add_argument('-e',type=str,dest="exon",help="STR BOOL. \"True\", for checking that a mate should overlap and exon to be counted as outgoing link of a window (if both mates are on the same window, at least one of them should overlap an exon). \"False\", if want to count all links.")
parser.add_argument('-w',type=int,dest="window",help="INT. Window size (nt) when splitting genome.")
parser.add_argument('-d',type=int,dest="dist",help="INT. Distance between mates to be consider not selfligating.")
parser.add_argument('-c',type=str,dest="cFile",help="TAB file. Chromosome size information")
parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file." )

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
awareFile = args.aFile
blindFile = args.bFile
checkExon = args.exon
windSize = args.window
dist = args.dist
chromFile = args.cFile
oFile = args.oFile



##################################
# General Functions

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

############################################
# Aware Links functions

def linkToIvAware(line,chromLength,windSize,dist):
# Create iv for DNA and RNA mates 
	line=line.split("\t")
	selfLig=True
	dna_iv = ivReg( HTSeq.GenomicInterval( line[0],int(line[1]),int(line[2]) ), chromLength, windSize )
	rna_iv = ivReg( HTSeq.GenomicInterval( line[9],int(line[10]),int(line[11]) ), chromLength, windSize )
	if line[16]=="exon" or line[16]=="utr3" or line[16]=="utr5": 
		coding_on=True
	else:
		coding_on=False
	if dna_iv.chrom!=rna_iv.chrom:
		selfLig=False
	elif abs(dna_iv.start-rna_iv.start)>dist:
		selfLig=False
	return [dna_iv,rna_iv,coding_on,selfLig]
	
def bed2awareLinks(File,chromLength,windSize,RNAonExon,dist):
# Convert DNA-RNA links from BED to GenomicArrayOfSets format
	links = {}
	
	for line in open(File,"r"):
		line = line.strip()	
		dna_iv, rna_iv, coding_on,selfLig = linkToIvAware(line,chromLength,windSize,dist)
		# Cosiders only RNA mates sitting on exons?
		if RNAonExon and (not coding_on): continue
		if selfLig: continue
			
		if rna_iv in links:
			if dna_iv in links[rna_iv]:
				links[rna_iv][dna_iv] += 1
			else:
				links[rna_iv][dna_iv] = 1
		else:
			links[rna_iv] = {}
			links[rna_iv][dna_iv] = 1

	return links


################################################
# Blind links functions

def lineToIvBlind(line,chromLength,windSize,dist):
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
	

def bed2blindLinks(bFile,chromLength,windSize,checkExon,dist):
# Convert DNA-RNA links from BED to GenomicArrayOfSets format

	links = {}
	
	for line in open(bFile,"r"):
		line = line.strip()	
		iv1, iv2, coding1,coding2,selfLig = lineToIvBlind(line,chromLength,windSize,dist)
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
awareLinks = bed2awareLinks(awareFile,chromLength,windSize,checkExon,dist) # rna_iv: set of dna_iv
blindLinks = bed2blindLinks(blindFile,chromLength,windSize,checkExon,dist) # rna_iv: set of dna_iv
print "Matching number of hits and peaks"
ratio={}
zeroc,allc=0,0
# Count hits on targets of region_iv
# Count histone peaks on region_iv
# Only counts peaks on regions that already have links counts
# This since most genomic regions having histone peaks aren't
# expected to have interaction with all aimers.				
for iv1 in awareLinks: # Iterate over aimer regions
	## Requieres that current aimer has at least 5 targets
	#if len(links[iv1])<5: continue

	# Count all blind-link targets
	allBlindCount=0
	if iv1 in blindLinks:
		for iv2 in blindLinks[iv1]:
			allBlindCount += blindLinks[iv1][iv2]

	# Count blind-links on aware-links targets
	specificBlindCount=0
	if iv1 in blindLinks:
		for iv2 in awareLinks[iv1]: # Iterate over targets of current aimer
			if iv2 in blindLinks[iv1]:
				specificBlindCount += blindLinks[iv1][iv2]
	
	allc+=1
	if allBlindCount==0:
		zeroc+=1
		continue
	ratio[iv1] = float(specificBlindCount) / float(allBlindCount)		

print "Printing results to output file"
out=open(oFile,'w')
for iv in ratio:
	print >>out, iv.chrom+"\t"+str(iv.start)+"\t"+str(iv.end)+"\t"+str(ratio[iv])
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
