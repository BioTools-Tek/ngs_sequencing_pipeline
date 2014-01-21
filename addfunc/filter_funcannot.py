#!/usr/bin/env python

import sys
from VCF_Operations import VCFOps

VERSION="0.7"
name = sys.argv[0].split('/')[-1]
ids="--ids="

def usage():
	print >> sys.stderr, '''ver %s
	
%s file.vcf %s<AL+DIL+COL+MUL+PRL+CCH+PCH> <[List of (Keyword|ID field|(ALL|N))>
                 or
%s --help  for examples''' % (VERSION, name, ids, name);
	exit(-1)

def help():
	print >> sys.stderr, '''
By default this script removes introns and intergenic regions where possible, but only if a keyword is given

e.g.  
%s file.vcf %s   SYNONYMOUS|MUL|1
    if no ids are given, then the script attempts to find the ids itself (95%% reliable)

%s file.vcf %sAL+DIL+COL+MUL+PRL+CCH+PCH SYNONYMOUS|MUL|1 Gly|PRL|1
   will extract just the lines where the the mutation field (identified by
   'MUL') have at least 1 SYNONYMOUS *OR* at least 1 protein Gly in the PRL
   field *OR* both.

%s file.vcf %sAL+DIL+COL+MUL+PRL+CCH+PCH SYNONYMOUS|MUL|ALL Gly|PRL|ALL
    same as above except that lines must have SYNONYMOUS in all elements
    of the field (ignoring '*') *OR* Gly in all elements of the PRL field
    *OR* both

%s file.vcf %sAL+DIL+COL+MUL+PRL+CCH+PCH SYNONYMOUS|MUL|1+Gly|PRL|ALL GGG|COL|2
    the '+' is an AND condition, interpreted as: Extract lines where there
    are 1 or more SYNONYMOUS mutations in the MUL fields *AND* all the elements
    in the PRL field are Gly *OR* 2 or more GGG in the COL field
 
Many AND and OR statements can be chained together using the '+' and ' ' seperators respectively
''' % (\
name, ids,\
name, ids,\
name, ids,\
name, ids
)
	exit(-1);

class Keyword:

	def __init__(self, key, field, num):
		self.keyword = key.strip()
		self.field = field.strip()
		self.foundtargets = 0                  # same as 'satisfied'
		
		num=num.strip()
		if (num=="ALL"):
			num = -100
		self.num = int(num);

class Condition:
	"Each condition is a collection of AND operators, and they interact with each other via ORring"
	"e.g. a condition with a single keyword will OR with another condition object"
	"a conditions with multiple keywords will AND the keywords, and OR with another condition object"
	def __init__(self, keywObjs):
		self.keyobjects = []
		self.satisfied = False
		
		objs = keywObjs.split('+')
		for ob in objs:
			osplit = ob.split('|')
			kw = Keyword(osplit[0],osplit[1],osplit[2])
			self.keyobjects.append(kw);


def parseArgs(args):
	
	lenarg = len(args)
	
	if lenarg==1:
		usage();
	
	if lenarg<4:
		if args[1]=="--help":
			help();
		usage();

	vcf_file = args[1];
	
	if args[2].startswith(ids):
		id_v=args[2].split(ids)[1];
	else:
		print >> sys.stderr, "Please give %s<something> as second argument, or " % ids
		exit(-1)
		
	condition_array = []
	for i in xrange(3,len(args)):
		c = Condition(args[i]);
		condition_array.append(c)
	
	return vcf_file, id_v.split('+'), condition_array



file, arg_ids, conditions = parseArgs(sys.argv);

prefix = file.split('/')[-1]
num_tabs = len(prefix)/8

vcf = VCFOps(file)
vcf.getFieldIndexes(arg_ids);

if arg_ids[0]=='':
	arg_ids = vcf.getSameLengthIds()
	print >> sys.stderr, "No Ids given, using:", arg_ids

print >> sys.stderr, prefix,

f = open(file,'r')

linecount=0
for line in f:
	linecount += 1
	
	line = line.splitlines()[0]
	if line[0]=='#':
		print line
		continue
		
	if linecount%11==0:
		print >> sys.stderr, '\r', num_tabs*'\t', '  #:', (100*linecount)/vcf.numlines,'%',
	
	tokens = line.split('\t')
	func_data = tokens[vcf.IFORMAT_INDEX].split(':')
	data = vcf.getElementsFromData(func_data, check=arg_ids)
	
	# Process arguments
	printline=False
	all_good_indexes=[]
	
	for cond in conditions:
		#OR at this level
		
		truth_list=[] #stores whether all keywords in condition have had their targets met

		for keyw in cond.keyobjects:
			#AND at this level
			index = vcf.field_map[keyw.field]
			list_of_data = data[keyw.field]
			

			max_freq = keyw.num
			len_data = len(list_of_data)
			
			good_indexes=[]  # stores indexes for that keyword
			num_empty_indexes=0   # stores frequency of unwanted indexes
			
			for k in xrange(len_data):
				a_key = list_of_data[k].strip()
				
				if a_key == "*":
					num_empty_indexes += 1
					continue
				
				if (a_key == keyw.keyword):
					good_indexes.append(k)
				
#				if len(good_indexes) >= max_freq:
#					break  # target number of keyword reached
				
			num_good_res = len(good_indexes)
			
			
			keyw.foundtargets = 0
			
			#For 'ALL', need to check that it doesnt count '*'
			if max_freq==-100:	#if checking for ALL
				
				if (num_good_res > 0) and ((num_good_res + num_empty_indexes)==len_data):   # good + empties = max only if all other fields are not this keyword
					keyw.foundtargets = 1
				
			else:									# checking for greater or equal to max_freq
				if num_good_res >= max_freq:
					keyw.foundtargets=1

# note to self: How to debug
# ./filter_funcannot.py <file> --ids=AL+DIL+COL+PRL+MUL "SYNONYMOUS|MUL|ALL" 2>&1 >out | grep "\[ 0 \]" | grep SYNONYMOUS
#
# count the number of that keyword should appear per line (if it exists), then check to see the output file emulates this


# DEBUG each keyword
#			print >> sys.stderr, '\n', list_of_data,' == ', num_good_res, '+', num_empty_indexes, '=', len_data, '[', keyw.foundtargets, '] max=', max_freq
			truth_list.append(keyw.foundtargets);
			all_good_indexes += good_indexes;

		#check all keywords for the same condition are equal
		all_equal = ((reduce(lambda x,y: x+y, truth_list))==len(truth_list))
#		print >> sys.stderr, truth_list, all_equal
		cond.satisfied=all_equal;
		if cond.satisfied:  # as long as one condition is met, then print line
			
			line_data = []
			
			for vcf_ids in vcf.FORMAT_IDS_ORDER:
				try:
					current_list = data[vcf_ids]
				except KeyError:
					line_data.append( func_data[vcf.field_map[vcf_ids]] ) 
					continue
			
				new_list=[]
				for index in all_good_indexes:
					new_list.append(current_list[index])
				line_data.append(','.join(new_list))
			
			tokens[vcf.IFORMAT_INDEX] = ':'.join(line_data)
			printline = True
			break;
	
	if printline:
		print '\t'.join(tokens)

print >> sys.stderr, '\r', num_tabs*'\t', "  #: 100 %"
