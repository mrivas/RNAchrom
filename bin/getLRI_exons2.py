import argparse
import sys
import pysam
import HTSeq
import numpy
import re

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Counts the number of long-range interactions per gene, but only for RNA-ends sitting on top of exons, utr3, or utr5 regions. Both DNA and RNA are counted per gene')
parser.add_argument('-a',type=str,dest="aFile",help="BED file. Output of annotateBAM.py, where the left and right ends correspond to DNA and RNA, respectively")
parser.add_argument('-g',type=str,dest="gFile",help="GTF file. Annotation file with biotype information")
parser.add_argument('-d',type=int,dest="distance",help="INT. Minimum distance (in nucleotides) between RNA-DNA ends to be considered as a long range interaccion. Default = 4000 b",default=4000)
parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file. Default = \"output\"",default="output")

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
aFile = args.aFile
gtfFile = args.gFile
distance  = args.distance
oFile = args.oFile

#aFile="/data2/yu68/bharat-interaction/TwoStep1024-2014-3-14/split_partner/TwoStep1024_fragment_paired_align.txt"
#gtfFile="/home/rivasas2/tools/genomes/mouse/mm9/Mus_musculus.NCBIM37.67.gtf.gz"
#distance=5000
#oFile="test.txt"
##############################################################
# Functions

def getGenes(gtfFile):
	# Get overla genomic regions per each annotated feature in the annotation file

	geneID_IV = {}
	gtf_file = HTSeq.GFF_Reader( gtfFile )

	for feature in gtf_file:
		ID = feature.attr["gene_id"]
		if ID in geneID_IV:
			if feature.iv.start < geneID_IV[ID].start:
				geneID_IV[ID].start = feature.iv.start
			if feature.iv.end > geneID_IV[ID].end:
				geneID_IV[ID].end = feature.iv.end
		else:
			geneID_IV[ID] = feature.iv
	
	
	# Get ID per genomic region
	geneIV_ID = HTSeq.GenomicArrayOfSets("auto",stranded=False)
	for geneID in geneID_IV:
		iv = geneID_IV[geneID]
		geneIV_ID[ iv ] += geneID
		
	return geneID_IV,geneIV_ID

def getGenomicIntervals(line,distance):
	line=line.strip("\n")
	fields=line.split("\t")
	iv_dna= HTSeq.GenomicInterval(fields[0],int(fields[1]),int(fields[2]),".")
	iv_rna= HTSeq.GenomicInterval(fields[6],int(fields[7]),int(fields[8]),".")
	longRange = False

	# if RNA-end is not sitting on top of an exon, utr3, or utr5 region: skips
	if not ( fields[10].split("|")[-2] == "exon" ):
		return [iv_dna,iv_rna,longRange]

	# if distance ==0, consider all reads for the analysis=> longRange as True
	if distance == 0:
		longRange = True
	elif iv_dna.chrom != iv_rna.chrom:
		longRange = True
	elif min( abs(iv_dna.start - iv_rna.end), (iv_rna.start - iv_dna.end) ) > distance:
		longRange = True

	return [iv_dna,iv_rna,longRange]

def getLongRangeInteractions(aFile,distance,geneIV_ID):
	# Gets the counts of long range interactions for each gene

	countsDNA = {}
	countsRNA = {}
	for line in open(aFile,"r"):
		iv_dna, iv_rna, longRange = getGenomicIntervals(line,distance)
		
		if not longRange: continue

		####################################
		# Count DNA ends
		iset = None
		for iv,step_set in geneIV_ID[iv_dna].steps():
			if iset is None:
				iset = step_set.copy()
			else:
				iset.intersection_update( step_set )
		if len( iset ) == 1:
			key = list(iset)[0]
			if key in countsDNA:
				countsDNA[ key ] += 1
			else:
				countsDNA[ key ] = 1
		####################################
		# Count RNA ends
		iset = None
		for iv,step_set in geneIV_ID[iv_rna].steps():
			if iset is None:
				iset = step_set.copy()
			else:
				iset.intersection_update( step_set )
		if len( iset ) == 1:
			key = list(iset)[0]
			if key in countsRNA:
				countsRNA[ key ] += 1
			else:
				countsRNA[ key ] = 1

	return [countsDNA, countsRNA]

def getBioType(gtfFile):
	bioType={}
	
	gtf_file = HTSeq.GFF_Reader( gtfFile )

	for feature in gtf_file:
		ID = feature.attr["gene_id"]
		bioType[ ID ] = feature.attr["gene_biotype"]
	
	return bioType

####################################################
# Excecution

print "Getting genes"
geneID_IV,geneIV_ID = getGenes(gtfFile)
print "Counting long range interactions per gene"
countsDNA, countsRNA = getLongRangeInteractions(aFile,distance,geneIV_ID) 
print "Getting biotypes"
bioType = getBioType(gtfFile)
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
