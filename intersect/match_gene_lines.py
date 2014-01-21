#!/usr/bin/env python

import sys

def usage():
	print >> sys.stderr, '''
v0.4  Extracts lines from a VCF file that match the genes(|exons) specified in a genelist
	
	usage: %s infile.vcf  genelist.txt --geneid=AL [OPTIONS]

	--geneid     argument must be given
	
	OPTIONS:
	--exons      remove exons from gene names
	--isoforms   (REMOVED: it extracts only the gene names given in genelist)
	''' % sys.argv[0].split('/')[-1]
	exit(-1)
	

if len(sys.argv)<4:
	usage()
	

exons=False
#isoforms=False

for arg in sys.argv[4:]:
	if arg == "--exons":
		exons=True
#	elif arg == "--isoforms":
#		isoforms=True
	else:
		usage()

#1.  Handle genelist
try:
	genes = open(sys.argv[2],'r')
except IOError:
	print >> sys.stderr, "[Error] Cannot open genelist", sys.argv[2]
	exit(-1)
	
genelist = []
for f in genes:
	f = f.strip()
	if not f in genelist:
		genelist.append(f)
	else:
		print >> sys.stderr, "Duplicate:", f

genes.close()
#genelist.sort() # sorting not neccesary, in this case it might actually slow things down


#2. Handle VCF file
# Grab gene identifier ( default is 'AL')
g_id = sys.argv[3].split("--geneid=")[1]

try:
	vcf = open(sys.argv[1],'r')
except IOError:
	print >> sys.stderr, "[Error] Cannot open VCF file", sys.argv[1]
	exit(-1)

m_count = 0
for line in vcf:
	line = line.splitlines()[0].strip()
	
	if line.startswith('#'):
		print line
	else:
		tokens = line.split('\t')
		
		try:
			g_index = tokens[-2].split(':').index(g_id)
			g_list = tokens[-1].split(':')[g_index].split(',')
		except (IndexError,ValueError):
			print >> sys.stderr, "Cannot find %s for line:\n%s" % (g_id, tokens[-2])
			print >> sys.stderr, "Gene index, List: %d, %s" % (g_index, g_list)
			exit(-1)
		
		#Gene1|Exon1, Gene2|ExonD
		for gene_name in g_list:
			gene_name = gene_name.strip();

			if not exons:
				gene_name = gene_name.split('|')[0].strip()

			if gene_name.startswith("ADAM12"):
				print >> sys.stderr, "before:", gene_name

#			if not isoforms:
#				gene_name = gene_name.split("-ISO")[0].strip()
				
#				if gene_name.startswith("ADAM12"):
#					print >> sys.stderr, "after:", gene_name
					
				
			
#			print >> sys.stderr, gene_name
			
			if gene_name in genelist:
				m_count += 1
				print line
				break  # needs to only match one gene in the list to print the line

vcf.close()

name=sys.argv[1].split('/')[-1]
nchar_m=len(str(m_count))
nchar_n=len(name)
nspaces=40-(nchar_m+nchar_n)


print >> sys.stderr, name, nspaces*' ', " ---->  matched  #:",m_count
