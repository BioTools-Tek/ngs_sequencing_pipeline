#!/usr/bin/env python

import sys

al_arg="--geneid="
isoforms = False
exons = False

def usage():
	print >> sys.stderr, '''
Extracts gene names (and exons) from vcf files, sorts and uniqs them
	
	%s %sAL infile.vcf [OPTIONS]
	
	Options:
	--isoforms		includes isoforms as seperate genes
	--exons			includes exons as seperate genes
''' % (sys.argv[0], al_arg)
	exit(-1)

if len(sys.argv)<3:
	usage()

AL= sys.argv[1].split(al_arg)[1]
file = sys.argv[2]


for arg in sys.argv[3:]:
	if arg == "--isoforms":
		isoforms = True
	elif arg == "--exons":
		exons = True
	else:
		usage()

try:
	f = open(file, 'r')
except IOError:
	print >> sys.stderr, "[Error] Cannot open", file
	exit(-1)

genelist = []

for line in f:
	line = line.splitlines()[0].strip()
	
	if not(line.startswith('#')):
		tokens = line.split('\t')
		try:
			al_index = tokens[-2].split(':').index(AL)
			al_list = tokens[-1].split(':')[al_index].split(',')
		except IndexError:
			print >> sys.stderr, "Cannot find %s for line:\n%s" % (AL, line)
			print >> sys.stderr, "zyg index, osity: %d, %s" % (al_index, al_list)
			exit(-1)
		
		#Gene1|Exon1, Gene2|ExonD
		for gene_name in al_list:
			gene_name = gene_name.strip()
			
			if not exons:
				gene_name = gene_name.split('|')[0].strip()
			
			if not isoforms:
				gene_name = gene_name.split("-ISO")[0].strip()
			
			if not gene_name in genelist:
				genelist.append(gene_name)

f.close()
genelist.sort();

for g in genelist:
	print g

