import argparse
import sys
import pysam
import HTSeq
import numpy
import re

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Get gene intervals (including introns,exons,utrs, etc) and ouputs in BED format')
parser.add_argument('-g',type=str,dest="gFile",help="GTF file. Annotation file")
parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of BED output file.")

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
gFile = args.gFile
oFile = args.oFile

##############################################################
# Functions

################################################################
def getGenes(gtfFile):
	# Get overla genomic regions per each annotated feature in the annotation file
	
	gtf = HTSeq.GFF_Reader(gtfFile)
	geneID_IV = {}

	for feature0 in gtf:
		feature=feature0
		# Add "chr" prefix to Ensembl chromosome names
		if not re.match("^chr",feature.iv.chrom):
			if feature.iv.chrom == "MT": feature.iv.chrom = "M"
			feature.iv.chrom = "chr"+feature.iv.chrom
		
		ID = feature.attr["gene_id"]+"|"+feature.attr["gene_name"]+"|"+feature.attr["gene_biotype"]
		if ID in geneID_IV:
			if feature.iv.start < geneID_IV[ID].start:
				geneID_IV[ID].start = feature.iv.start
			if feature.iv.end > geneID_IV[ID].end:
				geneID_IV[ID].end = feature.iv.end
		else:
			geneID_IV[ID] = feature.iv

	return geneID_IV



####################################################
# Excecution

print "Getting genes"
genes = getGenes(gFile)
print "Saving results to: "+oFile
out = open(oFile,'w')
for key in genes:
	iv = genes[key]
	output=[iv.chrom,iv.start,iv.end,key,".",iv.strand]
	print >>out, "\t".join(map(str,output))
			
out.close()
