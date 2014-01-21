#!/usr/bin/env python

# This is heavily based off filter_funcannot.py but I could not get the
# filtering to work without keywords, so seperate scripts

import sys
from time import sleep
from VCF_Operations import VCFOps

VERSION="0.2"
name = sys.argv[0].split('/')[-1]
ids="--ids="
noisof="--no-isoforms"
nospl="--no-splice"

def usage():
	print >> sys.stderr, '''ver %s

	usage: %s file.vcf %s<AL+DIL+COL+MUL+PRL+CCH+PCH> [%s] [%s]

To autofind ids, simply give '%s' and the script will attempt to find them.
This script removes introns by default, but can also remove isoforms and splice sites
''' % (VERSION, name, ids, noisof, nospl, ids);
	exit(-1)


def parseArgs(args):
	if len(args)<3:
		usage();

	noisoforms=False
	nosplice=False
	
	vcf_file = args[1].strip();
	
	if args[2].startswith(ids):
		id_v=args[2].split(ids)[1];
	else:
		print >> sys.stderr, "Please give %s<something> as second argument, or " % ids
		exit(-1)
		
	for i in xrange(3,len(args)):
		if args[i].startswith(noisof):
			noisoforms=True;
		elif args[i].startswith(nospl):
			nosplice=True;

	return vcf_file, id_v.split('+'), noisoforms, nosplice




file, arg_ids, noisoforms, nosplice = parseArgs(sys.argv);

prefix = file.split('/')[-1]
num_tabs = len(prefix)/8

linecount= 0
vcf = VCFOps(file)
vcf.getFieldIndexes(arg_ids);

if arg_ids[0]=='':
	arg_ids = vcf.getSameLengthIds()
	print >> sys.stderr, "No Ids given, using:", arg_ids


f = open(file,'r')
print >> sys.stderr, prefix,

for line in f:
	linecount += 1
	line = line.splitlines()[0]
	
	if line[0]=='#':
		print line
		continue
	
	tokens=line.split('\t')
	
	func_data = tokens[vcf.IFORMAT_INDEX].split(':')
	data_map = vcf.getElementsFromData(func_data, check=arg_ids);
	
	if linecount%11==0:
		print >> sys.stderr, '\r', num_tabs*'\t', '  #:', (100*linecount)/vcf.numlines,'%',
	
	all_bad_indexes=[]
	length_data = len(data_map[arg_ids[0]])
	
	
	for k in xrange(length_data):
		bad_indexes = []
		
		#Have length of data, now check each id for each
		for id_name in arg_ids:
			element = data_map[id_name][k];
			
			if element == "*":
				bad_indexes.append(k)
				continue
			
			if element.find("Intron")!=-1:
				bad_indexes.append(k)
				continue
			
			if noisoforms:
				if element.find("-ISOF")!=-1:
					bad_indexes.append(k)
					continue
			
			if nosplice:
				if element.find("Splice")!=-1:
					bad_indexes.append(k)
					continue
		
		all_bad_indexes += bad_indexes
	all_bad_indexes = list(set(all_bad_indexes)) #uniq
	
	if len(all_bad_indexes)==length_data:
		continue
	
	# Data gathered, now reassemble line
	# Reassembles in the order of the arg IDs, regardless of how it is in the VCF file	data = []
	data = []
	
	for vcf_ids in vcf.FORMAT_IDS_ORDER:
		try:
			current_list = data_map[vcf_ids]
		except KeyError:
			data.append( func_data[vcf.field_map[vcf_ids]] ) 
			continue
			
		new_list= []
		for index in xrange(length_data):
			if index not in all_bad_indexes:
				new_list.append(current_list[index])
		data.append ( ','.join(new_list) )

	tokens[vcf.IFORMAT_INDEX] = ":".join(data)
	
	print '\t'.join(tokens)

print >> sys.stderr, '\r', num_tabs*'\t', "  #: 100 %"
