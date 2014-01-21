#!/usr/bin/env python

from sys import stderr as cerr, argv

TYPE_id = "TYP"
SNV_id = "SNV"
IND_id = "IND"

HFORMAT = "##FORMAT="

def usage():
	print >> cerr, '''
Annotates the VCF file to specify which lines are indels and which are snvs denoted as "%s" and "%s" respectively

	usage: %s infile.vcf (snv|indel) > outfile.vcf
	''' % ("snv", "indel", argv[0])
	exit(-1)

if len(argv)<3:
	usage()

if not(argv[2]=="snv" or argv[2]=="indel"):
	print >> cerr, "Error:  Please specify either snv or indel", argv[2]
	usage()

vcf=argv[1]
IDS = SNV_id if (argv[2]=="snv") else IND_id

#Header implementation taken from bamzygo
def handleHeaders(file):
	
	HEADER_ID = HFORMAT+"<ID="+TYPE_id+','
	HEADER_ID_FULL = HEADER_ID+'Number=.,Type=String,Description="Marks line as single nucleotide variant(%s) or insertion/deletion(%s)">' % (SNV_id, IND_id)
	
	try:
		f = open(file,'r')
	except IOError:
		print >> cerr, "Could not open", file
		exit(-1)


	found_id = False
	found_format = -1

	for line in f:
		line = line.splitlines()[0].strip()
		
		if not line.startswith('#'):
			break
				
		#Look for ##
		if line.startswith(HFORMAT):
			found_format =0; #found
			if line.startswith(HEADER_ID):
				found_id = True;
			
		else:
			if (found_format==0):
				if not(found_id):
					print HEADER_ID_FULL
				found_format=-2
				
		print line
						
	f.close()


def appendToAll(file,id):
	try:
		f = open(file,'r')
	except IOError:
		print >> cerr, "Could not open", file
		exit(-1)

	for line in f:
		if line.startswith('#'):
			continue
		
		line = line.splitlines()[0].strip()
		tokens = line.split('\t')
		
		exists=False
		
		try:
			index1 = tokens[-2].split(':').index(TYPE_id)
		except ValueError:
			index1 = -1
			
		try:
			
			if index1==-1:	#Can't find ID, make one, index is now last
				tokens[-2] += ':'+TYPE_id
				#index = -1, already points to last index
			else:
				exists=True

			if not exists:
				tokens[-1] += ':'+id
			else:			#CAN find id, replace
				splitokens = tokens[-1].strip().split(':')
				splitokens[index1] = id
				tokens[-1] = reduce(lambda x,y: x+':'+y, splitokens)
		
		except IndexError:
			print >> cerr, "Problem:", exists, splitokens, index1
			exit(-1)
		
		print reduce(lambda x,y: x+'\t'+y, tokens)
	f.close()

handleHeaders(vcf)
appendToAll(vcf, IDS)
