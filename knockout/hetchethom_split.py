#!/usr/bin/env python

import sys
import collections

# This script seperates lines into three different data files:
# Homozygous, Single Heterozygous, and Compound Heterozygous

gene_id = "--geneid="
zygo_id = "--zygid="
typ_id = "--typid="
depth_q = "--depth="
ex_id = "--exons"
is_id = "--isoforms"
ignore_q = "--ignore-heteroq"
printlines_q = "--printlines"


#constants
delim="::::"
SNV="SNV"
IND="IND"
HET="HETEROZY"
HOM="HOMOZY"
HETQ="HETEROZY?"

def usage():
	print >> sys.stderr,'''
Splits lines of a VCF file into Homozygous, Heterozygous, and Compound Heterozygous genes for SNVs
	
	usage: %s <vcf-like file> <output_prefix> <[OPTIONS]>
	
Required:
%sAL	Identifier for the gene ID
%sZYG	Identifier for zygosity
%sTYP	Identifier for Type (for SNP and IND)

Optional:
[%s]	Print lines horizontally corresponding to each gene in the genelist
[%s]	Ignore 'HETEROZY?' lines. By default these are treated as 'HETEROZY'
[%s]	Use exon specific seperation of genes
[%s]	Use isoform specific seperation of genes

''' % (sys.argv[0].split('/')[-1], gene_id, zygo_id, typ_id, printlines_q, ignore_q, ex_id, is_id)
	exit(-1)

########################## Handle Args #################################
def parseargs():
	if len(sys.argv)<5:
		usage()

	################ Required ##########################################
	vcf=sys.argv[1]
	prefix=sys.argv[2]

	################ Optional ##########################################

	g_id = ""
	z_id = ""
	t_id= ""
	
	exons=False
	isoforms=False
	ignore_hetq=False
	printlines=False

	#args
	for arg in sys.argv[3:]:
		if arg.startswith(ex_id):
			exons=True
		elif arg.startswith(is_id):
			isoforms=True
		elif arg.startswith(ignore_q):
			ignore_hetq=True
		elif arg.startswith(printlines_q):
			printlines=True
		elif arg.startswith(gene_id):
			g_id=arg.split(gene_id)[1]
		elif arg.startswith(zygo_id):
			z_id=arg.split(zygo_id)[1]
		elif arg.startswith(typ_id):
			t_id=arg.split(typ_id)[1]

	######## Handle Errors #############################################

	if g_id == "" or z_id == "" or t_id == "":
		print >> sys.stderr, "Must specify zyg, gene, and type ids"
		print >> sys.stderr, g_id, z_id, t_id
		exit(-1)

	return vcf, prefix, g_id, z_id, t_id, exons, isoforms, ignore_hetq, printlines

########################################################################

########################################################################


def getHeaders(vcf):
	# Find headers (if any) close file, reopen
	headers=[]
	try:
		vfile = open(vcf,'r')
	except IOError:
		print >> sys.stderr, "Error:", vcf
		exit(-1)

	tmp = vfile.readline().splitlines()[0];
	headers.append(tmp);
	while (tmp.find("#CHROM")==-1):
		headers.append(tmp)
		tmp= vfile.readline().splitlines()[0];
	headers.append(tmp)
	vfile.close()
	return headers


def populateMap(vcf,g_id,z_id, t_id, exons, isoforms):
	gene_dict={}
	#[gene][zygosity] =  [line1,line27]  #num lines gives freq

	try:
		file = open(vcf,'r')
	except IOError:
		print >> sys.stderr, "Cannot open", vcf

	for line in file:
		line = line.splitlines()[0].strip()

		if not line.startswith('#'):
			tokens = line.split('\t')
			
			try:
				tokes2 = tokens[-2].split(':')

				g_index = tokes2.index(g_id)
				z_index = tokes2.index(z_id)
				t_index = tokes2.index(t_id)
				
			except IndexError:
				print >> sys.stderr, "temp"
				exit(-1)
				
			try:
				tokes1 = tokens[-1].split(':')
				
				g_list = tokes1[g_index].strip()
				zygosity = tokes1[z_index].strip().split(',')[-1]
				typ = tokes1[t_index].strip()
				
				if not(typ==SNV or typ==IND):
					print >> "Cannot parse line:", line
				
				
			except (IndexError,ValueError):
				print >> sys.stderr, "Cannot find %s for line:\n%s" % (g_id, tokens[-2])
				print >> sys.stderr, "Gene index, List: %d, %s" % (g_index, g_list)
				exit(-1)
			
			#Gene1|Exon1, Gene2|ExonD
			glist = g_list.split(',')
			for gene_name in glist:
				gene_name = gene_name.strip();
				
				# Don't count Introns or other unwanted 'genes'
				tmp = gene_name.split('|')
				try:
					if not(tmp[1].startswith("Exon") or tmp[1].startswith("Splice")):
						continue
				except IndexError:
					pass

				if not exons:
					gene_name = tmp[0].strip()

				if not isoforms:
					gene_name = gene_name.split("-ISO")[0].strip()

#				if gene_name == "XPC-ISOF2"

				#Map stuff
				if not gene_name in gene_dict:
					zygo_dict = {}
					zygo_dict[zygosity] = [line]

					type_map = {}
					type_map[typ] = zygo_dict

					gene_dict[gene_name] = type_map
					#                      __[HET]
					# Gene1 ____ [SNV] ___|__[HOM]
					#         |           |__[HET?]
					#         |__[IND] __ same as above

				elif not typ in gene_dict[gene_name]:
					zygo_dict = {}
					zygo_dict[zygosity] = [line]
					gene_dict[gene_name][typ] = zygo_dict
					
				elif not zygosity in gene_dict[gene_name][typ]:
					gene_dict[gene_name][typ][zygosity]= [line]
				
				elif not line in gene_dict[gene_name][typ][zygosity]:
					gene_dict[gene_name][typ][zygosity].append(line)
#
# Otherwise duplicate lines are printed, since the same gene name can occur more than once on the same line if we are ignoring isoforms, etc
#				else:
#					print >> sys.stderr, "DUPLICATE!"
#					exit(-1);

	return gene_dict


def printOut(prefix, headers, gene_dict, ignore_hetq, printlines):
	
	def mydict():
		return collections.defaultdict(mydict)
	
	try:
		snp="_snp"
		indel="_indel"
		
		__HET="het"
		__HOM="hom"
		__CHET="chet"
		__GENES="genes"
		__VCF="vcf"
		
		# Initialise file map
		file_map= {}
		for k1 in [SNV,IND,(SNV+IND)]:
			file_map[k1] = {}
			for k2 in [__HET,__HOM,__CHET]:
				if (k1==(SNV+IND) and not k2==__CHET):
					continue
				file_map[k1][k2] = {}
				for k3 in [__GENES,__VCF]:
					file_map[k1][k2][k3] = open(prefix+"_"+k1+"."+k2+"."+k3, 'w')


		#print headers
		for k1,v1 in file_map.iteritems():
			for k2,v2 in v1.iteritems():
				for k3,v3 in v2.iteritems():
					if k3==__VCF:
						print >> v3, '\n'.join(headers)

	except IOError:
		print >> sys.stderr, "Could not open output files for prefix=", prefix
		exit(-1)

	#                      __[HET]
	# Gene1 ____ [SNV] ___|__[HOM]
	#         |           |__[HET?]
	#         |__[IND] __ same as above
	#
	# gene_dict[gene_name][typ][zygosity].append(line)
	
	if not ignore_hetq:
		print "HETQ added to HET"


	print "GENE\tTYPE\t%s\t%s\t%s" % (HET,HOM,HETQ)
	print ""
	
	for gene in gene_dict.keys():
		type_map = gene_dict[gene]
		
		hetsnv_number=-1 #Stores het numbers across type
		hetind_number=-1
		hetsnv_lines = []
		hetind_lines = []
		
		
		for typ in type_map.keys():
			zyg_map = type_map[typ]
			
			print gene,'\t', typ,'\t',
			
			num_het = num_hom = num_hetq = 0
			
			######## HET ##################
			try:
				hetq_lines = zyg_map[HETQ]   # Find hetq lines first...
				num_chet = len(hetq_lines)
			except KeyError:
				pass
			
			try:
				het_lines = zyg_map[HET]
				num_het = len(het_lines)
				
				if not ignore_q:             #... so we can add them here if needed
					het_lines += hetq_lines
					num_het += num_chet
				
				if num_het == 1:
					print >> file_map[typ][__HET][__GENES], gene, ("\t%d\t%s" % (num_het, delim.join(het_lines))) if printlines else ""
					print >> file_map[typ][__HET][__VCF], '\n'.join(het_lines)
				elif num_het >= 2:
					print >> file_map[typ][__CHET][__GENES], gene, ("\t%d\t%s" % (num_het, delim.join(het_lines))) if printlines else ""
					print >> file_map[typ][__CHET][__VCF], '\n'.join(het_lines)
				else:
					print >> sys.stderr, "No het"
					
				# Store variables outside loop
				if typ == SNV:
					hetsnv_number = num_het
					hetsnv_array = het_lines
				elif typ == IND:
					hetind_number = num_het
					hetind_array = het_lines
				else:
					print >> sys.stderr, "NOTHING HERE"


			except KeyError:
				pass
			
			###### HOM #########################################
			try:
				hom_lines = zyg_map[HOM]
				num_hom = len(hom_lines)
				
				if num_hom >=1:
					print >> file_map[typ][__HOM][__GENES], gene, ("\t%d\t%s" % (num_hom, delim.join(hom_lines))) if printlines else ""
					print >> file_map[typ][__HOM][__VCF], '\n'.join(hom_lines)
				else:
					print >> sys.stderr, "No hom"
				
			except KeyError:
				pass
			
			print num_het,'\t',num_hom,'\t',num_hetq
		
		####### HET SNV and HET IND ##############
		# Find overlap between SNV and Indel het's for the same gene
		if hetsnv_number >=1 and hetind_number >=1:
			array = hetsnv_array + hetind_array
			array = sorted(set(array))
			number = hetsnv_number + hetind_number
			
			print >> file_map[SNV+IND][__CHET][__GENES], gene, ("\t%d\t%s" % (number, delim.join(array))) if printlines else ""	
			print >> file_map[SNV+IND][__CHET][__VCF], '\n'.join(array)



########################################################################
########################## M A I N #####################################

vcf, prefix, g_id, z_id, t_id, exons, isoforms, ignore_hetq, printlines = parseargs()

headers = getHeaders(vcf)
gene_dict = populateMap(vcf, g_id, z_id, t_id, exons, isoforms)

if gene_dict=={}:
	print >> sys.stderr, "Gene dictionary empty"
	exit(-1)

#print gene_dict['T'].keys()
#print gene_dict['T']['SNV'].keys()
#print gene_dict['T']['SNV']['HETEROZY']
#gene_dict[gene_name][typ][zygosity].append(line)

printOut(prefix, headers, gene_dict, ignore_hetq, printlines)
