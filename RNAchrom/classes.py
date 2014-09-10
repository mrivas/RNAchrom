class dotDict(dict):
	def __getattr__(self, attr):
		return self.get(attr, None)
	__setattr__= dict.__setitem__
	__delattr__= dict.__delitem__
########################################
class annotReader:
	
	def __init__(self, annotFile):
		self.annotFile = open(annotFile,'r')
		self.mate1, self.mate2 = dotDict,dotDict
		self.lineNumber = -1 # dummy initial value to be used on __getitem__'s "if" stament 

	def __getitem__(self,index):
		line = self.annotFile.readline().strip().split("\t")
		# Refresh annotFile for reuse. 
		if self.annotFile.tell()==self.lineNumber:
			self.annotFile.seek(0)
		self.lineNumber   = self.annotFile.tell()
		self.mate1.iv      = HTSeq.GenomicInterval(line[0],int(line[1]),int(line[2]),line[3] )
		self.mate1.genes   = line[4].split('|')[0].split(',')
		self.mate1.exon    = line[4].split('|')[1]
		self.mate1.repeat  = line[4].split('|')[2]
		self.mate1.sjs     = eval( line[5] )
		self.mate1.dist    = min( abs( int(line[2])-int(line[7]) ), abs( int(line[1])-int(line[8]) ) )
		self.mate1.readName = line[12]
		self.mate2.iv     = HTSeq.GenomicInterval(line[6],int(line[7]),int(line[8]),line[9] )
		self.mate2.genes  = line[10].split('|')[0].split(',')
		self.mate2.exon   = line[10].split('|')[1]
		self.mate2.repeat = line[10].split('|')[2]
		self.mate2.sjs    = eval( line[11] )
		self.mate2.dist   = min( abs( int(line[2])-int(line[7]) ), abs( int(line[1])-int(line[8]) ) )
		self.mate2.readName = line[12]
		return [self.mate1,self.mate2]	
		
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
for mate1,mate2 in awareMates:
	if not( len(mate2.genes)==1 and mate2.genes!=['.'] and mate2.exon=="exon" and mate1.sjs==[-1] and mate1.dist>dist ):
		continue 
	interactions[ mate2.genes[0]] = mates1.iv
	awareLinks[ mate1.readName  ] = [mate2.genes[0], mate1.iv ]
# Build blindLinks and add inferred aware-links to awareLinks
blindMates=RNAchrom.annotReader(blindAnnot)
blindLinks=RNAchrom.links()
for mate1,mate2 in blindMates:
	if not( len(mate2.genes)==1 and mate2.genes!=['.'] and mate2.exon=="exon" and mate1.sjs==[-1] and mate1.dist>dist ):
		continue 
	blindLinks[mate1.readName ] = [mate2.genes[0], mate1.iv  ]
	blindLinks[mate1.readName ] = [mate1.genes[0], mate2.iv ]

# Count hits on targets
for aimerGene, targetIV in interactions:
	awareCount,awareSD = awareLinks.count( aimerGene, targetIV )
	blindCount,blindSD = blindLinks.count( aimerGene, targetIV )
		
	print >>out, aimerGene, geneIV[aimerGene],codingSize[aimerGene],targetIV,awareCount,awareSD,blindCount,blindSD



