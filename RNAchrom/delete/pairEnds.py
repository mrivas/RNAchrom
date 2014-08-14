import argparse
import sys
import HTSeq

###########################################################
# Get command line arguments

parser = argparse.ArgumentParser(description='Output the positions of reads that are simultaneously mapped in files 1 and 2.')
parser.add_argument('-f1',type=str,dest="f1",help="STR. File 1 name.")
parser.add_argument('-f2',type=str,dest="f2",help="STR. File 2 name.")
parser.add_argument('-o',type=str,dest="oFile",help="STR. Name of output file. Default = \"output\"",default="output")

if len(sys.argv) == 1:
	print >> sys.stderr,parser.print_help()
	exit(0)

args = parser.parse_args()
f1 = args.f1
f2 = args.f2
oFile = args.oFile

##############################################################
# Execution

bam1 = HTSeq.BAM_Reader(f1)
bam2 = HTSeq.BAM_Reader(f2)

pair1={}
for almnt in bam1:
	if not almnt.aligned or almnt.not_primary_alignment: continue
	pair1[almnt.read.name] = almnt.iv

out=open(oFile,'w')

for almnt in bam2:
	if not almnt.aligned or almnt.not_primary_alignment: continue
	if almnt.read.name in pair1:
		iv1=pair1[almnt.read.name]
		print >>out, iv1.chrom+"\t"+str(iv1.start)+"\t"+str(iv1.end)+"\t"+iv1.strand+"\t"+almnt.read.name+"\t"+almnt.iv.chrom+"\t"+str(almnt.iv.start)+"\t"+str(almnt.iv.end)+"\t"+almnt.iv.strand

out.close()

