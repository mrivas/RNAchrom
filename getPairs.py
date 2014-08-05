import argparse
import sys
import HTSeq
import re

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Creates WashU interaction files from a paired-end BAM file.')
parser.add_argument('-b',type=str,dest="bFile",help="BAM file. Aligned paired-end file")
parser.add_argument('-o',type=str,dest="oFile",help="Output file.")

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
bamName = args.bFile
oFile = args.oFile

##############################################################
# Outputs the both mates in a single line of a BED file
bamFile = HTSeq.BAM_Reader(bamName)
out=open(oFile,"w")

first,second = None,None
count = 0

for almnt in bamFile:
	if not( almnt.aligned and almnt.paired_end and almnt.mate_aligned):
		continue
	if first == None:
		first = almnt
	else:
		second = almnt
		# Check that first and second have same read name
		if first.read.name != second.read.name:
			print "Pair names are different"
			print first.read.name, second.read.name
			first = None
			second = None
			continue
		
		if re.search("random",first.iv.chrom)  or re.search("random",second.iv.chrom):
			first = None
			second = None
			continue

		if first.iv.chrom == second.iv.chrom:
			if first.iv.start <= second.iv.start:
				strand1="+"
				strand2="-"
			else:
				strand1="-"
				strand2="+"
		else:
			strand1="."
			strand2="."
		
		count += 1
		output=[first.iv.chrom,first.iv.start,first.iv.end, second.iv.chrom+":"+str(second.iv.start)+"-"+str(second.iv.end)+",10",count,strand1]
		print >> out, '\t'.join(map(str,output))

		count += 1
		output=[ second.iv.chrom,second.iv.start,second.iv.end, first.iv.chrom+":"+str(first.iv.start)+"-"+str(first.iv.end)+",10",count,strand2]
		print >> out, '\t'.join(map(str,output))


		first = None
		second = None


out.close()
