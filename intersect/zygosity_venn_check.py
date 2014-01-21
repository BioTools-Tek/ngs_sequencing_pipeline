#!/usr/bin/env python

import sys
from time import sleep
from pprint import pprint

#Contstants
d_arg="--delim="
z_arg="--zygid="
i_arg="--ignore-lowhet"
a_arg="--agreement="

het="HETEROZY"
heter="HETEROZY?"
hom="HOMOZY"


delimterm="::::"
zygoterm="ZYG"

def usage():
	print >> sys.stderr, '''
Takes a venn output file (combined, overlapped) and filters out any lines that have inconsistent zygosities between individuals specified in the file
	
	%s <venn.output> <rejects.output> [options]"
	
where %s::::		is the delimeter used to seperate individuals on the same line
      %sZYG		is the keyword to find zygosity.
      %s 		ignores checking against HETEROZY? and assumes the line is correct
      %s100		is the percentage of agreement of zygosity across individuals
				the same position required for a 'good' line.
	
Defaults are displayed above, but must be given.
	''' % (sys.argv[0], d_arg, z_arg, i_arg, a_arg)
	exit(-1)


def parseArgs():
     lenarg = len(sys.argv)

     if lenarg < 2:
		 usage()

     #Defaults
     v_file = sys.argv[1]
     r_file = sys.argv[2]
     delim = delimterm
     zyg = zygoterm
     i_lowhet = False
     agree_lim=98


     index=lenarg-1
     while index > 0:
        arg = sys.argv[index]

        if arg.startswith('-'):
               
           if arg.startswith(i_arg):
               i_lowhet=True
               sys.argv.remove(arg)

           elif arg.startswith(d_arg):
               delim=sys.argv[index].split(d_arg)[1]
               sys.argv.remove(arg)
                  
           elif arg.startswith(z_arg):
               zyg=sys.argv[index].split(z_arg)[1]
               sys.argv.remove(arg)

           elif arg.startswith(a_arg):
               try:
                  agree_lim=sys.argv[index].split(a_arg)[1]
                  sys.argv.remove(arg)
                  agree_lim=int(agree_lim)
               except ValueError:
                  print >> sys.stderr, "[Error] Must give valid agreement limit"
                  exit(-1)

           else:
			   print >> sys.stderr, arg
			   usage()

        index -= 1

     return v_file, r_file, delim, zyg, i_lowhet, agree_lim



v_file, r_file, delim, zyg, i_lowhet, agree_lim = parseArgs()

try:
	f = open(v_file,'r')
	r = open(r_file,'w')
except IOError:
	print >> sys.stderr, "Cannot open files!"
	
for line in f:
	if len(line)<10:
		continue
	
	line = line.strip().splitlines()[0]
	
	#Print headers to both files
	if line.startswith(delim):
		print >> r, line
		print line
		continue
	elif line.startswith("Chr"):
		print >> r, line
		print line
		continue
	
	
	#At this stage only data remains
	
	zygomap = {}
	
	individuals = filter(lambda x: len(x)>0, line.split(delim)[1:])
	for data in individuals:
		tokens = data.strip().split('\t')
		try:
			zyg_index = tokens[-2].split(":").index(zyg)
			zygosity = tokens[-1].split(':')[zyg_index].split(',')[-1].strip()
		except IndexError:
			print >> sys.stderr, "Cannot find %s for line:\n%s [file=%s]" % (zyg, tokens[-2].split(':'), v_file)
			print >> sys.stderr, "zyg index, osity: %d, %s" % (zyg_index, zygosity)
			print >> sys.stderr, individuals
			exit(-1)

		if not(zygosity in zygomap):
			zygomap[zygosity] = 1
		else:
			zygomap[zygosity] += 1
		
	length = len(individuals)
	
	hom_perc = 0; het_perc = 0; hetq_perc = 0;
	# This needs to be done individually
	try:
		hom_perc = 100*zygomap[hom]/length
	except KeyError:
		hom_perc = 0
	try:
		het_perc = 100*zygomap[het]/length
	except KeyError:
		het_perc = 0
	try:
		hetq_perc = 100*zygomap[heter]/length
	except KeyError:
		hetq_perc = 0

	if i_lowhet:				#If ignoring "HET?", add percentages to other classes
		het_perc += hetq_perc;
		hom_perc += hetq_perc;
		
#	print >> sys.stderr, "###############"
#	print >> sys.stderr, line
#	pprint(zygomap, sys.stderr)
#	print >> sys.stderr, "hom=", hom_perc, "het=", het_perc, "hetq=", hetq_perc, "length=", length
	
	if ((het_perc >= agree_lim) or (hom_perc >= agree_lim)):
#		print >> sys.stderr, "HET/HOM"
		print line
	elif (hetq_perc >= agree_lim):
#		print >> sys.stderr, "HETQ"
		print line
	else:
#		print >> sys.stderr, "ELSE"
		print >> r, line
		
#	sleep(1)
