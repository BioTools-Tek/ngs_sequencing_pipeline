#!/usr/bin/env python
import sys

def usage():
	print >> sys.stderr, "Creates custom track that specifies every single variant in a VCF file -> BED.\nThe idea is to view the variants at each step of the pipeline"
	print >> sys.stderr, "\tusage: %s <VCF.file> <trackname>" % sys.argv[0].split('/')[-1]
	exit(-1);


if len(sys.argv)!=3:
	usage();

vcf_file=sys.argv[1]
track_name=sys.argv[2]


print "track name=\""+track_name+"\" description=\""+track_name+"\" visibility=1 itemRgb=\"On\""


for line in open(vcf_file,'r'):
	if line[0]=='#':
		continue
	
	tokens = line.split('\t');
	chrom = tokens[0]
	if chrom[0]!='c':
		chrom= "chr"+chrom;
	
	#Unrecognized by genome browser!
	if chrom=="chrMT":
		continue
	
	pos = int(tokens[1])
	name = tokens[2]
	
	print "%s\t%d\t%d\t%s" % (chrom, pos-1, pos , name)
