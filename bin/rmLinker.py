import numpy, itertools
from Bio.Seq import Seq

##########################################################
# Functions
def locate(queries,read,numMismatches):
	read =numpy.array( map(str,read)  )
	minMismatches=numMismatches+1 # dummy initial value
	matchLength=0   # dummy initial value
	output  = [ 0, len(queries.itervalues().next()), 0, len(read), 'none','none'] # startQ, endQ, startR, endR, strand, Type('match','trim','none')
	for orientation in queries: # Iterates over forward & reverse complemented queries
		query = queries[orientation]
		if orientation == "f" or orientation == "fc":
			strand = "+"
		elif orientation == "r" or orientation == "rc":
			strand = "-"
		for i in range( 0, len(read) ):
			startR = i
			endR   = min( len(read) , i+len(query) )
			startQ = 0
			endQ   = min( len(query), startQ	+(endR-startR) )
			# Allow only 2 mistmatches
			mismatches = sum(query[startQ:endQ]!=read[startR:endR])			   
			if (mismatches <= numMismatches) and (mismatches < minMismatches) and (endR-startR)>matchLength:
				minMismatches = mismatches
				matchLength   = endR - startR
				# Require at least 15 nt to call a match
				if matchLength>=15:
					output = [ startQ, endQ, startR, endR, strand, 'match' ]
				else: # If match happens in the last 15 nt assumed is just by chance and use it only for trimming
					output = [ startQ, endQ, startR, endR, strand, 'trim' ]				
	return output

def printMates(queries,mate1,mate2,numMismatches, detailsFile, mate1File, mate2File):

	read_id_f, read_seq_f, read_id2_f, read_qual_f = mate1
	read_id_f, read_seq_s, read_id2_s, read_qual_s = mate2
	startQ_f,endQ_f,startR_f,endR_f,strand_f, Type_f = locate(queries,read_seq_f,numMismatches)
	startQ_s,endQ_s,startR_s,endR_s,strand_s, Type_s = locate(queries,read_seq_s,numMismatches)
	# one or two ends with linker sequence at the beginning
	if startR_f<=15 or startR_s<=15:
		classification="ambiguous"
	# Coherent mates
	elif Type_f=="match" and Type_s=="match" and strand_f!=strand_s:
		classification="coherent_mates"
	# Un-coherent mates
	elif Type_f=="match" and Type_s=="match" and strand_f==strand_s:
		classification="uncoherent_mates"
	# only one matches
	elif (Type_f=="match" and Type_s!="match") or (Type_f!="match" and Type_s=="match"):
		classification="one_mate_linker"
	# none matches
	elif Type_f!="match" and Type_s!="match":   
		classification="no_linker"

	# Chop sequences and qualities
	read_seq_f, read_qual_f = read_seq_f[:startR_f], read_qual_f[:startR_f]
	read_seq_s, read_qual_s = read_seq_s[:startR_s], read_qual_s[:startR_s]   
	# Print details 
	details  = [ read_id_f ,startR_f,endR_f,strand_f,Type_f ]
	details +=		  [ startR_s,endR_s,strand_s,Type_s, classification ]
	print >> detailsFile, '\t'.join(map(str,details))
	# Print only non ambiguos fastq mates and attach classification to read_id
	if classification != "ambiguous":
		read_id = read_id_f + "__" + classification
		read_id2 = read_id2_f + "__" + classification
		print >> mate1File, '\n'.join(map(str,[read_id, read_seq_f, read_id2, read_qual_f]))
		print >> mate2File, '\n'.join(map(str,[read_id, read_seq_s, read_id2, read_qual_s]))


###########################################################################3

query="CTAGTAGCCCATGCAATGCGAGGA"
fastqFile1="/data2/sysbio/UCSD-sequencing/2014-3-26-Tri/CCGGRm_dupPE_RNA-RNA-03-26-2014_R1.fastq"
fastqFile2="/data2/sysbio/UCSD-sequencing/2014-3-26-Tri/CCGGRm_dupPE_RNA-RNA-03-26-2014_R2.fastq"
numMismatches=2

detailsFile=open("matesDetails.txt",'w')
mate1File  =open("mate1.fastq",'w')
mate2File  =open("mate2.fastq",'w')

# Excecution
# Define all possible orientations of the query
queries={} # stores forward and reverse complemented query
forward  = map(str,query)
forward_complement = map(str,Seq(query).complement())
reverse = map(str, str( query[::-1] ) )
reverse_complement = map(str, str( Seq(query).reverse_complement() ) )
queries["f"]=numpy.array( forward  )
queries["r"]=numpy.array( reverse )
queries["fc"]=numpy.array( forward_complement  )
queries["rc"]=numpy.array( reverse_complement )

# Loop over fastq reads
fastq1=open(fastqFile1,'r')
fastq2=open(fastqFile2,'r')
# Read mates
mate1=['dummy','dummy','dummy','dummy']
count=0
while mate1[0] != '':
	count += 1
	if count%100==0: print count
	mate1 = [fastq1.readline().strip() for i in range(4)]
	mate2 = [fastq2.readline().strip() for i in range(4)]
	# On mates, chop/trim linker sequence, classify, and print them
	printMates(queries,mate1,mate2,numMismatches,detailsFile,mate1File,mate2File)

fastq1.close()
fastq2.close()
detailsFile.close()
mate1File.close()
mate2File.close()

