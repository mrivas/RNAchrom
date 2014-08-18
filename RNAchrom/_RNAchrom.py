import HTSeq, sys, numpy, random, scipy.stats, pickle, ast

################################################################
def test():
	return "This is a test"

################################################################
def geneExonIV(gtfFile):
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
def repeats(rFile):
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
######################################################################
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
def chromLength(chromFile):
	""" Creates a DICT of chromosome lengths.

	:param chromFile: TXT file. TAB separated file where first and second columns correspond to chromosome name, and chromosome size (nt), respectively.
	:returns: DICT. Dictionary where keys and values are chrom names and chrom sizes, respectively.
	"""
	chromLength={}
	for line in open(chromFile,'r'):
		line=line.strip().split("\t")
		chromLength[line[0]] = int(line[1])
	return chromLength
################################################################
def ivReg(iv,chromLength,windSize):
	""" Regularization of intervals. The input iv is transformed to a canonical region. This to avoid partial overlapping, and to get well defined aimer and target windows.

	:param iv: HTSeq.GenomicInterval. Genomic interval to be regularized.
	:param chromLength: DICT. Dictionary of chromosome lengths.
	:param windSize: INT. Window size in which the genome is divided.
	:returns: HTSeq.GenomicInterval. Regularized genomic interval containing input iv. 
	"""
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
def bed2peaks(File,chromLength,windSize):
	"""Convert DNA-RNA links from BED to GenomicArrayOfSets format
	
	:param File: BED file. File name where ChIP-seq peaks are stored.
	:param chromLength: DICT. Dictionary of chromosome lengths.
	:param windSize: INT. Window size in which the genome is divided.
	:returns: DICT. Number of peaks per genomic window. 
	"""
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
def lineToIv(line,dist):
	""" Create iv for DNA and RNA mates

	:param line: STR. Line from annotatedBAM.py output.
	:param dist: INT. Distance between used as lower threshold to determine self-ligation.
	:returns: ARRAY. [HTSeq.GenomicInterval,HTSeq.GenomicInterval,BOOL,BOOL,BOOL]. Fist mate genomic interval, second mate genomic interval, 
	"""
	line=line.split("\t")
	coding1,coding2=False,False
	selfLig=True
	iv1 = HTSeq.GenomicInterval( line[0],int(line[1]),int(line[2]) )
	iv2 = HTSeq.GenomicInterval( line[6],int(line[7]),int(line[8]) )
	# check if mates overlap exons
	if line[4].split("|")[-2]=="exon": 
		coding1=True
	if line[10].split("|")[-2]=="exon": 
		coding2=True
	# Check if mates are selfligating
	if iv1.chrom!=iv2.chrom:
		selfLig=False
	elif min( abs(iv1.end-iv2.start), abs(iv1.start-iv2.end) )>dist:
		selfLig=False

	return [iv1,iv2,coding1,coding2,selfLig]
################################################################
def bed2links(bFile,chromLength,windSize,checkExon,dist,linkType):
	""" Convert DNA-RNA links from BED to GenomicArrayOfSets format
	
	:param bFile: TXT file. Annotated file of mates (output of annotateBAM.py).
	:param chromLength: DICT. Dictionary of chromosome lengths.
	:param windSize: INT. Window size in which the genome is divided.
	:param checkExon: BOOL. Whether to check if RNA-mates overlapp exons.
	:param dist: INT. Distance between used as lower threshold to determine self-ligation.
	:param linkType: STR. Type of mates: "aware", or "blind". 
	:returns: DICT[DICT]=INT. Dictionary of target per aimer. 	
	"""
	links = {}
	
	for line in open(bFile,"r"):
		line = line.strip()	
		iv1, iv2, coding1,coding2,selfLig = lineToIv(line,dist)
		iv1Reg=ivReg(iv1,chromLength,windSize)
		iv2Reg=ivReg(iv2,chromLength,windSize)
		# Ignore selfLigating links (mates closer to each other by less than 2k nt)
		if selfLig: continue

		if linkType=="aware":
			# Count only if iv2Reg overlap an exon
			if (not checkExon) or coding2:
				if iv2Reg in links:
					if iv1Reg in links[iv2Reg]:
						links[iv2Reg][iv1Reg] += 1
					else:
						links[iv2Reg][iv1Reg] = 1
				else:
					links[iv2Reg] = {}
					links[iv2Reg][iv1Reg] = 1
			
		else:
			if iv1Reg!=iv2Reg:
				# Count only if iv1Reg overlap an exon
				if (not checkExon) or coding1:
					if iv1Reg in links:
						if iv2Reg in links[iv1Reg]:
							links[iv1Reg][iv2Reg] += 1
						else:
							links[iv1Reg][iv2Reg] = 1
					else:
						links[iv1Reg] = {}
						links[iv1Reg][iv2Reg] = 1
				# Count only if iv2Reg overlap an exon
				if (not checkExon) or coding2:
					if iv2Reg in links:
						if iv1Reg in links[iv2Reg]:
							links[iv2Reg][iv1Reg] += 1
						else:
							links[iv2Reg][iv1Reg] = 1
					else:
						links[iv2Reg] = {}
						links[iv2Reg][iv1Reg] = 1
			else:
				# If iv1Reg==iv2Reg at least one end must overlap and exon
				if (not checkExon) or ( coding1 or coding2 ):
					if iv1Reg in links:
						if iv2Reg in links[iv1Reg]:
							links[iv1Reg][iv2Reg] += 1
						else:
							links[iv1Reg][iv2Reg] = 1
					else:
						links[iv1Reg] = {}
						links[iv1Reg][iv2Reg] = 1
	return links
##################################################################
def genes(gtfFile):
	""" Get overlap genomic regions per each annotated feature in the annotation file
	
	:param gtfFile: GTF file.
	:returns: [DICT,DICT]. Dictionaries of gene id (key) gene iv (value), and gene iv (key) and gene id (value), respectively. 
	"""
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
#############################################################
def counts(aFile,distance,geneIV_ID):
	""" Gets the counts of long range interactions for each gene
	
	:param aFile: Annotated file of aware links (output of annotateBAM.py).
	:oaran distance: INT. Distance between mates used a lower threshold to be call as not self-ligating.
	:param genIV_ID: DICT. Dictionary where keys and values are gene genomic intervals and gene ids, respectively.
	:returns: [DICT,DICT]. Array of dictionaries, corresponding to the count of DNA and RNA per gene id (key), respectively. 
	"""
	countsDNA = {}
	countsRNA = {}
	for line in open(aFile,"r"):
		iv_dna, iv_rna, coding1, coding2, selfLig = lineToIv(line,distance)
		
		if selfLig or (not coding2): continue
		####################################
		# Count DNA ends
		iset = None
		for iv,step_set in geneIV_ID[iv_dna].steps():
			if iset is None:
				iset = step_set.copy()
			else:
				iset.intersection_update( step_set )
		if len( iset ) >= 1:
			for key in iset:
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
		if len( iset ) >= 1:
			for key in iset:
				if key in countsRNA:
					countsRNA[ key ] += 1
				else:
					countsRNA[ key ] = 1

	return [countsDNA, countsRNA]
######################################################
def getBioType(gtfFile):
	"""Creates a dictionary of biotypes per gene.

	:param gtfFile: GTF file.
	:returns: DICT. Keys and values are gene id and biotypes, respectively
	"""
	bioType={}
	
	gtf_file = HTSeq.GFF_Reader( gtfFile )

	for feature in gtf_file:
		ID = feature.attr["gene_id"]
		bioType[ ID ] = feature.attr["gene_biotype"]
	
	return bioType

##########################################################
# detectLinks
##########################################################
def loadMates(annotFile,Type,dist):
	"""Creates a dictionary of mates. If input mates are blind, it infers aware-mates based on splicing information.

	:param annotFile: TXT file. Annotation file (output of annotateBAM.py).
	:param Type: STR. Type of links on annotFile. Could be: "aware" (second mate correspond to RNA-end), or "blind" (it's not known what mates correspond to the RNA-end).
	:param dist: INT. Distance (nt) between mates to be considered not self-ligating (mates closer than this value are discarded).
	:returns: DICT. Dictionary of mates, where the key is the geneID of the aimer, and the value, the interval of the target region. If Type=="aware", it output a single DICT of aware mates. If type=="blind", it output and array of two DICTs, the first are the inferred aware-mates, and the second the blind-mates.
	"""
	aware,inferred,blind={},{},{}
	for line in open(annotFile):
		line=line.strip().split('\t')
		# Skip self-ligating mates
		if min( abs(int(line[2])-int(line[7])), abs( int(line[1])-int(line[8]) )) > dist: continue
		readName = line[12]    # read name
		# If aware always first and second mates correspond to DNA and RNA, respectively 
		if Type=="aware":
			# Check consistency: DNA mates not spliced
			DNA_SJs = ast.literal_eval(line[5])
			spliced = False
			for DNA_SJ in DNA_SJs:
				if DNA_SJ >=0: spliced=True
			if spliced : continue
			# Check ambiguity: RNA end overlap one and only one gene
			genes=line[10].split('|')[0].split(',')
			if len(genes)>1 or genes[0]==".": continue
			aware[ readName ] = {} # initialized dict
			aware[ readName ][ 'geneName' ] = genes[0] # mate2 annotation
			aware[ readName ][ 'targetIv' ] = HTSeq.GenomicInterval(line[0],int(line[1]),int(line[2]) ) # mate1 iv = dna iv
		# Infer aware from blind-links
		elif Type=="blind":
			# Mates where one end has a known spliced junction, and the other is not spliced are inferred as aware-links
			Iv1_SJs = ast.literal_eval(line[5])
			Iv2_SJs = ast.literal_eval(line[11])
			spliced1,spliced2=False,False
			for Iv1_SJ in Iv1_SJs:
				if Iv1_SJ >=0: spliced1=True
			for Iv2_SJ in Iv2_SJs:
				if Iv2_SJ >=0: spliced2=True
			# Inferred aware-links
			if (spliced1 + spliced2)==1: 
			
				if spliced1:   # Iv1 is the RNA-end
					# Check ambiguity: RNA end overlap one and only one gene
					genes=line[4].split('|')[0].split(',')
					if len(genes)!=1 or genes[0]==".": continue
					inferred[ readName ] = {}
					inferred[ readName ][ 'geneName' ] = genes[0]
					inferred[ readName ][ 'targetIv' ] = HTSeq.GenomicInterval(line[6],int(line[7]),int(line[8]) ) # mate2 iv = dna iv
				elif spliced2: # Iv2 is the RNA-end
					# Check ambiguity: RNA end overlap one and only one gene
					genes=line[10].split('|')[0].split(',')
					if len(genes)!=1 or genes[0]==".": continue
					inferred[ readName ] = {}
					inferred[ readName ][ 'geneName' ] = genes[0]
					inferred[ readName ][ 'targetIv' ] = HTSeq.GenomicInterval(line[0],int(line[1]),int(line[2]) ) # mate1 iv = dna iv
			# Blind-links
			elif (spliced1 + spliced2)==0: # cases where both mates have sj are dicarded as ambiguous
				# Annotate iv1 as rna iv (aimer)
				genes=line[4].split('|')[0].split(',')
				if len(genes)!=1 or genes[0]==".": continue
				blind[ readName ] = {}
				blind[ readName ][ 'geneName' ] = genes[0]
				blind[ readName ][ 'targetIv' ] = HTSeq.GenomicInterval(line[6],int(line[7]),int(line[8]) ) # mate2 iv = dna iv
				# Annotate iv2 as rna iv (aimer)
				genes=line[10].split('|')[0].split(',')
				if len(genes)!=1 or genes[0]==".": continue
				blind[ readName ] = {}
				blind[ readName ][ 'geneName' ] = genes[0]
				blind[ readName ][ 'targetIv' ] = HTSeq.GenomicInterval(line[0],int(line[1]),int(line[2]) ) # mate1 iv = dna iv

	if Type=="aware":
		return aware
	elif Type=="blind":
		return [inferred,blind]

############################################################
def expandIv(iv,windSize):
	"""Expand genomic interval

	:param iv: HTSeq.GenomicInterval. Genomic interval to be expanded
	:param windSize: INT. Final size of the expanded iv
	:returns: HTSeq.GenomicInterval. Genomic interval centered at iv and with length windSize
	"""
	halfWind = round( float(windSize)/2 )
	iv.start = max(iv.start - halfWind,0)
	iv.end   = iv.end   + halfWind
	return iv
############################################################
def std(readName_ivs):
	"""Computes standard deviation of a list of genomic intervals

	:param iv_array: ARRAY[ HTSeq.GenomicInterval ]. Array of genomic intervals.
	:returns: INT. Standard deviation of input genomic intervals.
	"""
	mindpoints = []
	for readName_iv in readName_ivs:
		iv=readName_iv.split("_")[1]
		start=iv.split("-")[0]
		end=iv.split("-")[1]
		midpoints.append( numpy.mean([start,end]) )
	return numpy.std(midpoints)
############################################################
def buildLinks( awareMates,blindMates, windSize,genesID_IV, oFile ):
	"""Count hits of aware and blind-mates. It determine target areas based on aware-mates. 

	:param awareMates: DICT. Dictionary of aware mates.
	:param blindMates: DICT. Dictionary of blind mates.
	:param windSize: INT. Window size of expanded target mates.
	:param genes: DICT. Dictionary of genes' genomic intervals, where keys are gene ids and values HTSeq.GenomicIntervals.
	:param oFile: STR. Name of output TAB separated file, where results will be stored. The output format is: geneName geneid_iv target_iv awareCounts awareSD blindCounts blindSD.
	"""
	# Define targets using aware-links as seeds
	targetAreas = {}
	for readName in awareMates:
		geneName = awareMates[readName]['geneName']
		targetIvExp = expandIv( awareMates[readName]['targetIv'], windSize)

		# Determine targets' areas
		if not ( geneName in targetAreas ):
			targetAreas[ geneName ] = HTSeq.GenomicArrayOfSets("auto",stranded=False)
			targetAreas[ geneName ][ targetIvExp ] += "hit"
		else:
			targetAreas[ geneName ][ targetIvExp ] += "hit" 
	
	awareHits,blindHits= {},{}
	# Create dictionary of genomic array of sets for aware links
	for readName in awareMates:
		geneName = awareMates[readName]['geneName']
		targetIv = awareMates[readName]['targetIv']
		if not ( geneName in awareHits ):
			awareHits[ geneName ] = HTSeq.GenomicArrayOfSets("auto",stranded=False)
			awareHits[ geneName ][ targetIv ] += readName+"_"+str(targetIv.start)+"-"+str(targetIv.end)
		else:
			awareHits[ geneName ][ targetIv ] += readName+"_"+str(targetIv.start)+"-"+str(targetIv.end) 
	# Create dictionary of genomic array of sets for blind links
	for readName in blindMates:
		geneName = blindMates[readName]['geneName']
		targetIv = blindMates[readName]['targetIv']
		if not ( geneName in blindHits ):
			blindHits[ geneName ] = HTSeq.GenomicArrayOfSets("auto",stranded=False)
			blindHits[ geneName ][ targetIv ] += readName+"_"+str(targetIv.start)+"-"+str(targetIv.end)
		else:
			blindHits[ geneName ][ targetIv ] += readName+"_"+str(targetIv.start)+"-"+str(targetIv.end) 
	
	out = open(oFile,'w')
	# Count and print hits on target areas
	for geneName in targetAreas:
		for target_iv,target_val in targetAreas[geneName].steps():
			if len(target_val)==0: continue # skip empty (no seeded) areas
			# Count number of aware hits on target_iv
			iset = None
			for awareHits_iv,awareHits_val in awareHits[geneName][target_iv].steps():
				if iset is None: iset = awareHits_val.copy()
				else:            iset.itersection_update( awareHits_val )
			awareCounts = len( iset )
			awareSD    = std( iset )
			# Count number of blind hits on target_iv
			iset = None
			for blindHits_iv,blindHits_val in blindHits[geneName][target_iv].steps():
				if iset is None: iset = blindHits_val.copy()
				else:            iset.itersection_update( blindHits_val )
			blindCounts = len( iset )
			blindSD    = std( iset )
			# Print results
			geneID=geneName.split('-')[0]
			output = [geneName,genesID_IV[geneID],target_iv,awareCounts,awareSD,blindCounts,blindSD]
			print >>out, '\t'.join(map(str,output))
	out.close()



