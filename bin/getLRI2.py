import argparse
import sys
import pysam
import HTSeq
import numpy
import re

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Counts the number of long-range interactions per gene. Both DNA and RNA are counted per gene')
parser.add_argument('-a',type=str,dest="aFile",help="BED file. Output of Stitch-seq_Aligner, where the left and right ends correspond to DNA and RNA, respectively")
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
	
		# Add "chr" prefix to Ensembl chromosome names
		if not re.match("^chr",geneID_IV[ID].chrom):
			if geneID_IV[ID].chrom == "MT": geneID_IV[ID].chrom = "M"
			geneID_IV[ID].chrom = "chr"+geneID_IV[ID].chrom

	
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
	iv_rna= HTSeq.GenomicInterval(fields[9],int(fields[10]),int(fields[11]),".")
	longRange = False
	
	# if distance ==0, consider all reads for the analysis=> longRange as True
	if distance == 0:
		longRange = True
	elif iv_dna.chrom != iv_rna.chrom:
		longRange = True
	elif (iv_dna.start - iv_rna.end)>=distance or (iv_rna.start - iv_dna.end)>=distance:
		longRange = True

	return [iv_dna,iv_rna,longRange]

def getLongRangeInteractions(aFile,distance,geneIV_ID):
	# Gets the counts of long range interactions for each gene

	countsDNA = {}
	countsRNA = {}
	for line in open(aFile,"r"):
		iv_dna, iv_rna, longRange = getGenomicIntervals(line,distance)
		
		if longRange==False: continue

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
