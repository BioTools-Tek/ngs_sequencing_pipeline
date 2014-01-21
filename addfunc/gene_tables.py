#!/usr/bin/env python

import sys
from VCF_Operations import VCFOps

#args
g_arg="--geneid="
m_arg="--mutid="
c_arg="--codid="
p_arg="--protid="
z_arg="--zygid="
ne_arg="--no-exons"
fr_arg="--frequencies"
ty_arg="--types"
ex_arg="--kidney-expression"

def usage():
	print >> sys.stderr, '''
ver0.5
Takes in multiple functionally annotated VCF files and prints an HTML table showing the mutations that occured for each gene for each patient

	usage: %s %sAL %sMUL %sCCH %sPCH %sZYG [OPTIONS] [<VCF Files>]

OPTIONS:
%s	groups genes by gene name and not exons
%s	repeated variations are grouped and tallied, not duplicated
%s	groups transitions by types: (ins, del, snv) (mis,syn,non)
''' % (\
sys.argv[0].split('/')[-1], g_arg, m_arg, c_arg, p_arg, z_arg,\
ne_arg, fr_arg, ty_arg)
	exit(-1)


## Parse Args ##
if len(sys.argv)<7:
	usage()

g_id=sys.argv[1].split(g_arg)[1]
m_id=sys.argv[2].split(m_arg)[1]
c_id=sys.argv[3].split(c_arg)[1]
p_id=sys.argv[4].split(p_arg)[1]
z_id=sys.argv[5].split(z_arg)[1]

if ( "" in (g_id, m_id, c_id, p_id, z_id ) ):
	print >> sys.stderr, "Require Ids in the order given:\n"
	usage();

noexons=False
frequencies=False
types=False
files=[]

for arg in sys.argv[6:]:
	if arg.startswith(ne_arg):
		noexons=True
	elif arg.startswith(fr_arg):
		frequencies=True
	elif arg.startswith(ty_arg):
		types=True
	else:
		files.append(arg)

########################################################################
print >> sys.stderr, "Making map:",
########################################################################
def getType(codon_change):
	if (codon_change.find("ins")!=-1):
		return "ins"
	if (codon_change.find("del")!=-1):
		return "del"
	return "snv"

numtabs=2;
gene_map={}      # Populated across all files
                 # map -> Gene -> File -> List of codon and protein changes

arg_ids = [g_id, m_id, c_id, p_id]   #z_id omitted because it doesn't have the same number of elements

for file in files:
	
	file=file.strip()
	prefix=file.split('/')[-1].split('_')[0]

	vcf = VCFOps(file)
	vcf.getFieldIndexes(arg_ids+[z_id])
	
	f = open(file,'r')
	
	countlines=0

	for line in f:
		countlines += 1
		
		if line[0]=='#':
			continue
		
		func_data = line.split('\t')[vcf.IFORMAT_INDEX].split(':')
		data_map = vcf.getElementsFromData(func_data, check=arg_ids);
		
		zygosity = func_data[vcf.field_map[z_id]].split(',')[-1][0:3]
		
		for index in xrange(len(data_map[g_id])):
			
			gene = data_map[g_id][index]
			mutation = data_map[m_id][index]
			codon_change = data_map[c_id][index]
			protein_change = data_map[p_id][index]
			
			if not types:
				uniq = (codon_change, protein_change)
			else:
				uniq = ( getType(codon_change), mutation[0:3], zygosity )

			if noexons:
				gene = gene.split('|')[0]

			if gene in gene_map:
				if file in gene_map[gene]:
					gene_map[gene][file].append( uniq )
				else:
					gene_map[gene][file] = [ uniq ]
			else:
				gene_map[gene] = {}
				gene_map[gene][file] = [ uniq ]

		if (countlines%21==0):
			print >> sys.stderr, '\r', '\t'*numtabs, prefix,':', int(100*countlines/vcf.numlines), '%',
	
	print >> sys.stderr, '\r', '\t'*numtabs, "%s:100%%   " % prefix,
	f.close()
	numtabs += 1



########################################################################
print >> sys.stderr, "\nPrinting Table:",
########################################################################
genecount=0;                        # progress
numgenes = len(gene_map.keys())     # num genes

field_sep=','
element_sep=':'

###HEADERS
print 'Patient:',field_sep,
for file in files:
	print file.split('/')[-1],field_sep,
print 'All Same?\n\nGenes:'

for gene in gene_map:
	print gene,field_sep,
	
	all_same=[]  # stores build or format for each file across the gene and sees if they're equal
	
	for file in files:
		try:
			if frequencies:
				
				freq_table={}

				for uniq in gene_map[gene][file]:
					if uniq not in freq_table:
						freq_table[uniq] = 1
					else:
						freq_table[uniq] += 1
				
				formatted = [ "  ["+str(v)+"] "+reduce(lambda x,y: x+" "+y, key)+"  "+element_sep for key,v in freq_table.iteritems()]
				print reduce(lambda x,y: x+y, formatted)[2:-1], field_sep,
				
				if formatted not in all_same:
					all_same.append(formatted)

			else:
				build=""
				for key in gene_map[gene][file]:
					build += "  "+reduce(lambda x,y: x+" "+y, key)+"  "+element_sep
				build = build[2:-1]

				if build not in all_same:
					all_same.append(build)

				print build,field_sep,

		except KeyError:
			print "",field_sep, #file doesnt have that gene
	
	
	if len(all_same)==1: # no duplicates, same mutations across individuals (i.e. all same)
		print 'X',


	genecount += 1
	print ''
	print >> sys.stderr, '\r\t\t', int(100*genecount/numgenes), '%',
print >> sys.stderr, "\r\t\t100%  "

