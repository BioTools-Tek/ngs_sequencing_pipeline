#!/usr/bin/env python

import sys
from pprint import pprint

def usage():
	print >> sys.stderr, '''
	Takes in two or more gene lists and prints the overlap
	
	%s <[FILES]>
	''' % sys.argv[0]

if len(sys.argv)==1:
	usage()

files = sys.argv[1:]

uniq_map = {};

for file in files:
	try:
		f = open(file,'r')
	except IOError:
		print >> sys.stderr, "[Error] no such file", file
		exit(-1)
	for line in f:
		line = line.splitlines()[0].strip()
		if line not in uniq_map:
			uniq_map[line] = 1;
		else:
			uniq_map[line] += 1;
	
	f.close()
	
#print len(uniq_map)

num_files = len(files)
common_set = filter(lambda x: uniq_map[x]==num_files, uniq_map.keys())

#pprint(common_set)
genes = sorted(common_set)

for g in genes:
	print g
