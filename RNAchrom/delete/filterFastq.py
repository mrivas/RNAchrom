#!/usr/bin/env python
from RNAchrom import *

###########################################################
# Get command line arguments
def getArgs():
	parser = argparse.ArgumentParser(description='Filter out reads given a list of read IDs.')
	parser.add_argument('-l',type=str,dest="lFile",help="TXT file. List of fastq reads to be discarded. Column_1=readID; column_10=type of ID.")
	parser.add_argument('-f',type=str,dest="fFile",help="FASTQ file. Fastq file to be filter.")
	parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output FASTQ file.")

	if len(sys.argv) == 1:
		print >> sys.stderr,parser.print_help()
		exit(0)
	return parser
####################################################
def main():
	args = gerArgs().parse_args()
	listFile = args.lFile
	fastqFile = args.fFile
	oFile = args.oFile
	#listFile="/data2/yu68/bharat-interaction/TwoStep1024-2014-3-14/split_partner/fragment_align_detail.txt"
	#fastqFile="/data2/sysbio/UCSD-sequencing/2014-3-13-Bharat/TwoStep1024/Rm_dupPE_TwoStep1024_S2_R1.fastq"
	#oFile="test.fastq"

	print "Loading ID of reads to be eliminated"
	linker={}
	nLine=0
	for line in open(listFile,'r'):
		nLine+=1
		if nLine==1: continue
		ID=line.strip().split("\t")[0]
		linker[ "@"+ID ]=1

	print "Filtering out reads"
	out=open(oFile,'w')
	nLine=0
	count=4
	for line in open(fastqFile,'r'):
		line=line.strip()
		nLine+=1
		if nLine%4==1:
			ID=line.split(" ")[0]
			if not( ID in linker):
				print >>out, line
				count=1
		elif count<4:
			count += 1
			print >>out, line
	out.close()
####################################################
if __name__ == '__main__':
	main()
