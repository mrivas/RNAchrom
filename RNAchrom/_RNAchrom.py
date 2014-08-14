import HTSeq
import sys
import pysam
import numpy
import re
import random
import scipy.stats
import pickle

################################################################
def test():
	return "This is a test"

################################################################
# ANNOTATEBAM
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

################################################################
# corrWind
################################################################
def getchromLength(chromFile):
	chromLength={}
	for line in open(chromFile,'r'):
		line=line.strip().split("\t")
		chromLength[line[0]] = int(line[1])
	return chromLength
################################################################
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
################################################################
def bed2Peaks(File,chromLength,windSize):
# Convert DNA-RNA links from BED to GenomicArrayOfSets format
	coverage = {}
	
	for line in open(File,"r"):
		line = line.strip().split("\t")	
		iv = ivReg( HTSeq.GenomicInterval( line[0],int(line[1]),int(line[2]) ), chromLength, windSize )
		if iv in coverage:
			coverage[iv] += 1
		else:
			coverage[iv] = 1

	return coverage
################################################################
def lineToIv(line,chromLength,windSize,dist):
# Create iv for DNA and RNA mates 
	line=line.split("\t")
	coding1,coding2=False,False
	selfLig=True
	iv1 = ivReg( HTSeq.GenomicInterval( line[0],int(line[1]),int(line[2]) ), chromLength, windSize )
	iv2 = ivReg( HTSeq.GenomicInterval( line[6],int(line[7]),int(line[8]) ), chromLength, windSize )
	# check if mates overlap exons
	if line[4].split("|")[-2]=="exon": 
		coding1=True
	if line[10].split("|")[-2]=="exon": 
		coding2=True
	# Check if mates are selfligating
	if iv1.chrom!=iv2.chrom:
		selfLig=False
	elif abs(iv1.start-iv2.start)>dist:
		selfLig=False

	return [iv1,iv2,coding1,coding2,selfLig]
################################################################
def bed2Links(bFile,chromLength,windSize,checkExon,dist,linkType):
# Convert DNA-RNA links from BED to GenomicArrayOfSets format

	links = {}
	
	for line in open(bFile,"r"):
		line = line.strip()	
		iv1, iv2, coding1,coding2,selfLig = lineToIv(line,chromLength,windSize,dist)
		# Ignore links not present on known chromosomes
		if type(iv1)!=HTSeq._HTSeq.GenomicInterval or type(iv2)!=HTSeq._HTSeq.GenomicInterval: continue
		# Ignore selfLigating links (mates closer to each other by less than 2k nt)
		if selfLig: continue

		if linkType=="aware":
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
#############################################################################
# longRangeInterac.py
#############################################################################
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
###################################################################	
# specificity
###################################################################	
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

################################################
# Blind links functions

def lineToIv(line,chromLength,windSize,dist):
# Create iv for DNA and RNA mates 
	line=line.split("\t")
	coding1,coding2=False,False
	selfLig=True
	iv1 = ivReg( HTSeq.GenomicInterval( line[0],int(line[1]),int(line[2]) ), chromLength, windSize )
	iv2 = ivReg( HTSeq.GenomicInterval( line[6],int(line[7]),int(line[8]) ), chromLength, windSize )
	# check if mates overlap exons
	if line[4].split("|")[-2]=="exon": 
		coding1=True
	if line[10].split("|")[-2]=="exon": 
		coding2=True
	# Check if mates are selfligating
	if iv1.chrom!=iv2.chrom:
		selfLig=False
	elif min( abs(iv1.end-iv2.start), abs(iv1.start-iv2.end) ) > dist:
		selfLig=False

	return [iv1,iv2,coding1,coding2,selfLig]
	
def bed2links(bFile,chromLength,windSize,checkExon,dist,linkType):
# Convert DNA-RNA links from BED to GenomicArrayOfSets format

	countL=0
	links = {}
	
	for line in open(bFile,"r"):
		line = line.strip()
		iv1, iv2, coding1,coding2,selfLig = lineToIv(line,chromLength,windSize,dist)
		# Ignore selfLigating links (mates closer to each other by less than 2k nt)
		if selfLig: continue
		
		# For aware links, iv2=rna_end, iv1=dna_end
		if linkType=="aware":
			if (not checkExon) or coding2:
				if iv2 in links:
					if iv1 in links[iv2]:
						links[iv2][iv1] += 1
					else:
						links[iv2][iv1] = 1
				else:
					links[iv2] = {}
					links[iv2][iv1] = 1
		# For blind links
		else:
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

