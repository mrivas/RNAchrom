import argparse
import sys
import pysam
import HTSeq
import numpy
import re

def getArgs():
	# Get command line arguments
	parser = argparse.ArgumentParser(description='Annotate both mates of paired-end alignments')
	parser.add_argument('-b',type=str,dest="bFile",help="BAM file(s). A single aligment file contains all paired-end mates, or a comma separated list of 2 files, one per each set of mates (useful when mate origin is important).")
	parser.add_argument('-g',type=str,dest="gFile",help="GTF file. Annotation file with biotype information")
	parser.add_argument('-r',type=str,dest="rFile",help="BED file. Annotation file with repeats information")
	parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file. Default = \"output\"",default="output")

	if len(sys.argv) == 1:
		print >> sys.stderr,parser.print_help()
		exit(0)
	return parser

###########################################################
def main():
	args = getArgs().parse_args()
	bFile = args.bFile
	gFile = args.gFile
	rFile = args.rFile
	oFile = args.oFile
	# Excecution #####################################################
	print "Getting repeats"
	repeat = getRepeats(rFile)
	print "Getting genes,biotype, and exons"
	gene,exon = getGenes(gFile)
	print "Saving results to: "+oFile
	out = open(oFile,'w')
	first = None
	firstMate={}
	for i in range( len(bFile.split(",") ) ):
		bam=HTSeq.BAM_Reader( bFile.split(",")[i] )
		
		for almnt in bam:
			###########################################################
			# If one file with both mates
			if len(bFile.split(","))==1:
				# Discard not aligned reads, and multihits
				if (not almnt.aligned) or almnt.optional_field("NH")>1: continue
				if first == None:    
					first = almnt
				else: 
					second = almnt
					# Check that first and second reads have same read name
					if first.read.name != second.read.name:
						first=second
						continue
					annotFirst=annotate(first,gene,exon,repeat)
					annotSecond=annotate(second,gene,exon,repeat)
					
					if first.pe_which=="first": # Always print first in pair at the beginning
						output=formatOutput(first,second,annotFirst,annotSecond)
					else:
						output=formatOutput(second,first,annotSecond,annotFirst)
					print >>out,output  
					first = None
			###########################################################
			# If two files, each with one mate
			elif len(bFile.split(","))>1:
				# Discard not aligned reads, and multihits
				if (not almnt.aligned) or almnt.optional_field("NH")>1: continue
				annot=annotate(almnt,gene,exon,repeat)
				
				if i==0:
					firstMate[almnt.read.name]=[almnt]
					firstMate[almnt.read.name].append( annot )
				elif i>0 and (almnt.read.name in firstMate):
					output=formatOutput(firstMate[almnt.read.name][0],almnt,firstMate[almnt.read.name][1],annot)
					print >>out, output

	out.close()

################################################################
def getGenes(gtfFile):
	"""Stores genomic intervals of genes (including introns) and exons in dictionaries.
	
	:param gtfFile: GTF file.
	:returns: [gene,exon] GenomicArrayOfSets, gene[genomic_iv]=set([geneID1,geneID2,...]), GenomicArrayOfSets, exon[genomic_iv]="exon"
	"""
	gtf_file = HTSeq.GFF_Reader( gtfFile )
	geneID_IV = {}
	exon    = HTSeq.GenomicArrayOfSets("auto",stranded=False)
	for feature0 in gtf_file:
		feature=feature0
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
	"""Stores genomic intervals of repeats as HTSeq.GenomicArrayOfSets.
	
	:param rFile: BED file with the name of repeats on 4th column.
	:returns: HTSeq.GenomicArrayOfSets[genomic_iv]=set([repeatID])
	"""
	repeat = HTSeq.GenomicArrayOfSets("auto",stranded=False)
	for line in open(rFile):
		line=line.strip().split("\t")
		iv=HTSeq.GenomicInterval(line[0],int(line[1]),int(line[2]))
		repeat[iv]+=line[3]
	return repeat

#######################################################################		
def annotate(almnt,genes,exon,repeat):
	"""Annotate an alignment
	
	:param almnt: Alignment, HTSeq.GenomicInterval
	:param genes: Genes, gene[genomic_iv]=set([geneID1,geneID2,...])
	:param exon: Exons, exon[genomic_iv]=set(["exon"])
	:param repeat: Repeats regions, HTSeq.GenomicArrayOfSets[genomic_iv]=set([repeatID])
	:returns: Annotation string, formed as: "[gene_annot1,gene_annot2,...|.]|[exon,intron|.]|[repeat_annot|.]". Dot (.) denotes absence of annotation.
	"""
	geneAnnot=set()
	exonAnnot=set()
	repeatAnnot=set()

	for cigop in almnt.cigar:
		if cigop.type!="M": continue
		# Gene level annotationG
		for iv,geneIDs in genes[cigop.ref_iv].steps():
			geneAnnot |= geneIDs
		# Exon annotation
		for iv,exons in exon[cigop.ref_iv].steps():
			exonAnnot |= exons
		# Repeat annotation
		for iv,repeats in repeat[cigop.ref_iv].steps():
			repeatAnnot |= repeats

	if len(geneAnnot)==0: geneAnnot="."
	else:                 geneAnnot=",".join(map(str,geneAnnot))
	if len(exonAnnot)==0: exonAnnot="."
	else:                 exonAnnot="exon"
	if geneAnnot!="." and exonAnnot==".": exonAnnot="intron"
	if len(repeatAnnot)==0: repeatAnnot="."
	else:                   repeatAnnot=",".join(map(str,repeatAnnot))
	
	return geneAnnot+"|"+exonAnnot+"|"+repeatAnnot

####################################################
def formatOutput(first,second,annotFirst,annotSecond):
	""" Takes alignments' information and creates an string fit to be printed.
	
	:param first: Genomic interval, HTSeq.GenomicInterval
	:param second: Genomic interval, HTSeq.GenomicInterval
	:param annotFirst: Annotation of *first* 
	:param annotSecond: Annotation of *second*
	:returns: String of annotations, "first.chrom, first.start, first.end, first.strand, first.annot, first.jMField, second.chrom,..."
	"""
	output=[first.iv.chrom,str(first.iv.start),str(first.iv.end),first.iv.strand]
	output=output+[annotFirst,first.optional_field("jM")]
	output=output+[second.iv.chrom,str(second.iv.start),str(second.iv.end),second.iv.strand]
	output=output+[annotSecond,second.optional_field("jM")]
	output=output+[first.read.name]
	output="\t".join(map(str,output))
	return output

####################################################
if __name__ == '__main__':
	main()
