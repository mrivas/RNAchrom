import argparse
import sys
import pysam
import HTSeq
import numpy
import re

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Annotate both mates of paired-end alignments')
parser.add_argument('-b',type=str,dest="bFile",help="BAM file(s). A single aligment file contains all paired-end mates, or a comma separated list of 2 files, one per each set of mates (useful when mate origin is important).")
parser.add_argument('-g',type=str,dest="gFile",help="GTF file. Annotation file with biotype information")
parser.add_argument('-r',type=str,dest="rFile",help="BED file. Annotation file with repeats information")
parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file. Default = \"output\"",default="output")

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
bFile = args.bFile
gFile = args.gFile
rFile = args.rFile
oFile = args.oFile


#bFile="/data2/yu68/bharat-interaction/TwoStep1024-2014-3-14/RDI_TwoStep1024_total_sort.bam"
#gFile="/home/rivasas2/tools/genomes/mouse/mm9/Mus_musculus.NCBIM37.67.gtf.gz"
#rFile="/home/rivasas2/tools/genomes/mouse/mm9/repeatMasker.UIUC.bed"
#oFile="test.bed"
##############################################################
# Functions

################################################################
def getGenes(gtfFile):
	# Get overla genomic regions per each annotated feature in the annotation file

	gtf_file = HTSeq.GFF_Reader( gtfFile )
	geneID_IV = {}
	exon    = HTSeq.GenomicArrayOfSets("auto",stranded=False)

	for feature0 in gtf_file:
		feature=feature0
		# Add "chr" prefix to Ensembl chromosome names
		if not re.match("^chr",feature.iv.chrom):
			if feature.iv.chrom == "MT": feature.iv.chrom = "M"
			feature.iv.chrom = "chr"+feature.iv.chrom
		
		ID = feature.attr["gene_id"]+"-"+feature.attr["gene_biotype"]
		if ID in geneID_IV:
			if feature.iv.start < geneID_IV[ID].start:
				geneID_IV[ID].start = feature.iv.start
			if feature.iv.end > geneID_IV[ID].end:
				geneID_IV[ID].end = feature.iv.end
		else:
			geneID_IV[ID] = feature.iv

		if feature.type=="exon":
			exon[feature.iv] += "exon"	
		
	
	# Get ID per genomic region
	geneIV_ID = HTSeq.GenomicArrayOfSets("auto",stranded=False)
	biotypeIV_ID = HTSeq.GenomicArrayOfSets("auto",stranded=False)
	for geneID in geneID_IV:
		iv = geneID_IV[geneID]
		geneIV_ID[ iv ] += geneID
		
	return [geneIV_ID,exon]

#######################################################################
def getRepeats(rFile):
	repeat = HTSeq.GenomicArrayOfSets("auto",stranded=False)
	for line in open(rFile):
		line=line.strip().split("\t")
		iv=HTSeq.GenomicInterval(line[0],int(line[1]),int(line[2]))
		repeat[iv]+=line[3]
	return repeat

#######################################################################		
def annotate(almnt,genes,exon,repeat):
	geneAnnot=[]
	exonAnnot=[]
	repeatAnnot=[]
	
	# Gene level annotation
	for iv,geneIDs in genes[almnt.iv].steps():
		if len(geneIDs)==0: continue
		for geneID in geneIDs:
			geneAnnot.append(geneID)
	if len(geneAnnot)==0: geneAnnot="."
	else:                 geneAnnot=",".join(map(str,geneAnnot))
	# Exon annotation
	for iv,geneIDs in exon[almnt.iv].steps():
		if len(geneIDs)==0: continue
		for geneID in geneIDs:
			exonAnnot.append(geneID)
	if len(exonAnnot)==0: exonAnnot="."
	else:                 exonAnnot="exon"
	if geneAnnot!="." and exonAnnot==".": exonAnnot="intron"
	# Repeat annotation
	for iv,geneIDs in repeat[almnt.iv].steps():
		if len(geneIDs)==0: continue
		for geneID in geneIDs:
			repeatAnnot.append(geneID)
	if len(repeatAnnot)==0: repeatAnnot="."
	else:                   repeatAnnot=",".join(map(str,repeatAnnot))

	return geneAnnot+"|"+exonAnnot+"|"+repeatAnnot

####################################################
def formatOutput(first,second,annotFirst,annotSecond):
	output=[first.iv.chrom,str(first.iv.start),str(first.iv.end),first.iv.strand]
	output=output+[annotFirst]
	output=output+[second.iv.chrom,str(second.iv.start),str(second.iv.end),second.iv.strand]
	output=output+[annotSecond]
	output=output+[first.read.name]
	output="\t".join(map(str,output))
	return output


####################################################
# Excecution

print "Getting repeats"
repeat = getRepeats(rFile)
print "Getting genes,biotype, and exons"
gene,exon = getGenes(gFile)
print "Saving results to: "+oFile
out = open(oFile,'w')
###########################################################
first = None
firstMate={}
for i in range( len(bFile.split(",") ) ):
	bam=HTSeq.BAM_Reader( bFile.split(",")[i] )
	
	for almnt in bam:
		# Discard not aligned reads, and multihits
		if (not almnt.aligned) or almnt.optional_field("NH")>1: continue
	
		###########################################################
		# If one file with both mates
		if len(bFile.split(","))==1:
			if first == None:    
				first = almnt
			else: 
				second = almnt
				# Check that first and second reads have same read name
				if first.read.name != second.read.name:
					#print "Pair names are different"
					#print first.read.name, second.read.name
					first=second
					continue
				#else:
					#print first.read.name, second.read.name

				annotFirst=annotate(first,gene,exon,repeat)
				annotSecond=annotate(second,gene,exon,repeat)

				output=formatOutput(first,second,annotFirst,annotSecond)
				print >>out,output  
				
				first = None

		###########################################################
		# If two files, each with one mate
		elif len(bFile.split(","))>1:
			annot=annotate(almnt,gene,exon,repeat)
			
			if i==0:
				firstMate[almnt.read.name]=[almnt]
				firstMate[almnt.read.name].append( annot )
			elif i>0 and (almnt.read.name in firstMate):
				output=formatOutput(firstMate[almnt.read.name][0],almnt,firstMate[almnt.read.name][1],annot)
				print >>out, output

out.close()
