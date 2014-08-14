annotateBAM
==========
*1 getGenes
getRepeats
annotate
formatOutput

corrWind
=======
*3 getchromLength
*4 ivReg
bed2Peaks
*2 lineToIv
*5 bed2Links

longRangeInterc
===============
*1 getGenes
*2 getGenomicIntervals
getLongRangeInteractions
getBioType

specificity
===========
*3 getchromLength
*4 ivReg
*2 linToIv
*5 bed2Links

def getGenes(gtfFile): #on longRangeInterc
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
		iv = geneID_	VIV[geneID]
		geneIV_ID[ iv ] += geneID

	return geneID_IV,geneIV_ID

