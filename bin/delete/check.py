from RNAchrom import *

###########################################################
# Get command line arguments

def getArgs():
	parser = argparse.ArgumentParser(description='Save to output file only the unmapped reads that have a linker sequence on them.')
	parser.add_argument('-u',type=str,dest="uFile",help="SAM file. List of unmapped reads who may have the linker sequence on them.")
	parser.add_argument('-f',type=str,dest="fFile",help="FASTQ file. Fiel containing either RNA or DNA ends.")
	parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output FASTQ file.")

	if len(sys.argv) == 1:
		print >> sys.stderr,parser.print_help()
		exit(0)
	return parser

####################################################
def main():
	args = getArgs().parse_args()
	unmapFile = args.uFile
	fastqFile = args.fFile
	oFile = args.oFile


	unmap={}
	nLine=0
	for line in open(unmapFile,'r'):
		nLine+=1
		if nLine%4!=1: continue
		ID=line.strip().split("\t")[0].split('/')[0]
		unmap[ ID ]=1

	out=open(oFile,'w')
	nLine=0
	count=4
	for line in open(fastqFile,'r'):
		line=line.strip()
		nLine+=1
		if nLine%4==1:
			ID=line.split(" ")[0]
			if ID in unmap:
				print >>out, line
				count=1
		elif count<4:
			count += 1
			print >>out, line
	out.close()

####################################################
if __name__ == '__main__':
	main()
