class annotReader:
	
	def __init__(self, annotFile):
		self.annotFile = open(annotFile,'r')
		self.lineNumber = -1 # dummy initial value to be used on __getitem__'s "if" stament 

	def __getitem__(self,index):
		line = self.annotFile.readline().strip().split("\t")
		# Refresh annotFile for reuse. 
		if self.annotFile.tell()==self.lineNumber:
			self.annotFile.seek(0)
		self.lineNumber   = self.annotFile.tell()
		self.firstIV      = HTSeq.GenomicInterval(line[0],int(line[1]),int(line[2]),line[3] )
		self.firstGenes   = line[4].split('|')[0].split(',')
		self.firstExon    = line[4].split('|')[1]
		self.firstRepeat  = line[4].split('|')[2]
		self.firstSJs     = eval( line[5] )
		self.secondIV     = HTSeq.GenomicInterval(line[6],int(line[7]),int(line[8]),line[9] )
		self.secondGenes  = line[10].split('|')[0].split(',')
		self.secondExon   = line[10].split('|')[1]
		self.secondRepeat = line[10].split('|')[2]
		self.secondSJs    = eval( line[11] )
		self.dist         = min( abs( int(line[2])-int(line[7]) ), abs( int(line[1])-int(line[8]) ) )
		self.readName     = line[12]
		return self	
		
#########################################	
class interactions():
	def __init__(self,windSize):
		self.link = {}
		self.windSize = windSize

	def __setitem__(self,gene,iv):
		if not( gene in self.link ):
			self.link[gene] = HTSeq.GenomicArrayOfSets("auto",stranded=False)
		midpoint=numpy.mean([iv.start,iv.end])
		iv.start = max(0,midpoint - round(float(self.windSize)/2) )
		iv.end   =       midpoint + round(float(self.windSize)/2)
		self.link[gene][iv] += 1
	
	def storeAsArray(self):
		self.genes= []
		self.targets = []
		for gene in self.link:
			for iv,val in self.link[gene].steps():
				if len(val)>0:
					self.genes.append(gene)
					self.targets.append(iv)

	def __getitem__(self,i):
		if i==0: self.storeAsArray()
		if i<len(self.genes):
			return [self.genes[i],self.targets[i]]
		else:
			raise StopIteration

#########################################	
class links:
	def __init__(self):
		self.links = HTSeq.GenomicArrayOfSets("auto",stranded=False)
	def __setitem__(self,readName,tup):
		gene,iv = tup
		coord=str(iv.start)+"-"+str(iv.end)
		self.links[iv] += readName+"_"+coord+"_"+gene # We add readName as suffix to avoid discard different hist from the same gene
	def filterReadsByGene(self,reads,aimerGene):
		x=set([])
		if len(reads)==0:
			return x
		for read in reads:
			readName,coord,gene = read.split("_")
			if gene == aimerGene:
				x.add( read )
		return x
	def getSD(self,uset):
		"""Computes standard deviation of a list of genomic intervals

		:param iv_array: ARRAY[ HTSeq.GenomicInterval ]. Array of genomic intervals.
		:returns: INT. Standard deviation of input genomic intervals.
		"""
		midpoints = []
		for read in uset:
			start,end=read.split("_")[1].split("-")
			start,end=int(start),int(end)
			midpoints.append( numpy.mean([start,end]) )
		return numpy.std(midpoints)
	def getCount(self,aimerGene,targetIV):
		uset = set([])	# union set
		for iv,reads in self.links[targetIV].steps():
			uset |= self.filterReadsByGene(reads,aimerGene)
		count = len(uset)
		if count!=0:
			sd     = self.getSD(uset)
		else:
			sd     = "nan"
		return [count,sd]
#########################################	

# Build target regions and awareHits
awareMates=RNAchrom.annotReader(awareAnnot)
awareLinks=RNAchrom.links()
interactions=RNAchrom.interactions(windSize)
for mates in awareMates:
	if not( len(mates.secondGenes)==1 and mates.secondGenes!="." and mates.secondExon=="exon" and mates.firstSJs==-1 and mates.dist>dist ):
		continue 
	interactions[ mates.secondGenes[0] ] = mates.firstIV
	awareLinks[   mates.readName       ] = [mates.secondGenes[0], mates.firstIV ]
# Build blindLinks and add inferred aware-links to awareLinks
blindMates=RNAchrom.annotReader(blindAnnot)
blindLinks=RNAchrom.links()
for mates in blindMates:
	if not( len(mates.secondGenes)==1 and mates.secondGenes!="." and mates.secondExon=="exon" and mates.firstSJs==-1 and mates.dist>dist ):
		continue 
	blindLinks[mates.readName ] = [mates.secondGenes[0], mates.firstIV  ]
	blindLinks[mates.readName ] = [mates.firstGenes[0] , mates.secondIV ]

# Count hits on targets
for aimerGene, targetIV in interactions:
	awareCount,awareSD = awareLinks.count( aimerGene, targetIV )
	blindCount,blindSD = blindLinks.count( aimerGene, targetIV )
		
	print >>out, aimerGene, geneIV[aimerGene],codingSize[aimerGene],targetIV,awareCount,awareSD,blindCount,blindSD



