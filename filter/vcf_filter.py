#!/usr/bin/env python

import sys, re

name = sys.argv[0]
AL = "AL"

def usage():
	print >> sys.stderr, '''
%s input.vcf [OPTIONS]      v0.2 (base: protofilter)

 OPTIONS
 -f <effect> or --filter <effect1> [ <effectN> ]
     specifies optional keywords to filter.
     (e.g SYNONYMOUS, 3PRIME). Case sensitive regex.

 -u or --unfilter <effect1> [ <effectN> ]
     specifies an optional keyword to NOT include with
     results (e.g. UTR). Case sensitive regex.

  By default the filter matches ANY keyword specified.
  (e.g. A B, will match any line with A or B)
  To match ALL keywords use the '+' symbol to chain keywords
  together (e.g A+B, will match any line with both A and B)

 -d <Num> or --depth <Num>
     Only items with a read depth above the specified
     number are handled. Default is 0.
     
 -k  Keep novel

''' % name
	sys.exit(-1)


def parse_args():
	lenarg = len(sys.argv)

	#Populate args
	files=[]; ## More than one?

	#Default Args
	filt=unfilt=[]
	depth=0
	standard=False
	keepnovel=False
	
	index = lenarg-1

	while index > 0:
		arg=sys.argv[index]
	
		if arg[0]=="-":
			if arg=="-s" or arg=="--standardize" or arg=="--standardise":
				standard=True
				sys.argv.remove(arg)
				
			elif arg=="-k" or arg=="--keepnovel":
				keepnovel=True
				sys.argv.remove(arg);


			elif arg=="-f" or arg=="--filter":
				try:
					u_index = index+1
					while (u_index < lenarg) and (sys.argv[u_index][0] != "-"):
						u_index += 1
					filt=sys.argv[index:u_index]
				except IndexError:
					print >> sys.stderr, "Filter Error: Please give filter keyword"
					sys.exit(-1)

			elif arg=="-u" or arg=="--unfilter":
				try:
					u_index = index+1
					while (u_index < lenarg) and (sys.argv[u_index][0] != "-"):
						u_index += 1
					unfilt=sys.argv[index:u_index]
				except IndexError:
					print >> sys.stderr, "Unfilter Error: Please give unfilter keyword"
					sys.exit(-1)

			elif arg=="-d" or arg=="--depth":
				try:
					num = sys.argv[index+1]
					depth=int(num)
				except ValueError:
					print >> sys.stderr, "Depth Error: Please give a number"
					sys.exit(-1)
				except IndexError:
					print >> sys.stderr, "Depth Error: Please give a number"
					sys.exit(-1)
				sys.argv.remove(num)
				sys.argv.remove(arg)

			else:
				usage()

		index -=1

	#Clear already caught tokens
	sys.argv.remove(name)
	if len(filt)>0:
		for fect in filt:
			sys.argv.remove(fect) 
	if len(unfilt)>0:
		for fect in unfilt:
			sys.argv.remove(fect) 

	# assign remainder
	files = sys.argv
	usage() if len(files) == 0 else 1

	return files, filt[1:], unfilt[1:], depth, standard, keepnovel


def openFile(filen, flag='r'):
	try:
		f=open(filen,flag)
	except IOError:
		print >> sys.stderr, "File Error: Unable to open %s" % filen
	return f


def process(file, filter_name="", unfilter_name="", depth_lim=20, standard=False, keepnovel=False):

	def grabHeaderLine(wb):
		line= ""
		#Scan first 30 rows or so
		for line in wb:
			if line.find("CHROM")!=-1:
				break;
			print line.splitlines()[0];
		return line;

	questionable_lines = [] # Holds lines with questionable content
	variant_lines = []	#Holds entire matching line to be standardized later

	f = openFile(file)
	headers = grabHeaderLine(f).strip()
	if standard:
		variant_lines.append(headers)
	else:
		print headers


	for line in f:
		if (line[0]=='#') or (len(line) < 10):
			continue

		tokens1 = line.split('\t');

		CHROM = tokens1[0];
		POS =tokens1[1]
		ID = tokens1[2]
		REF= tokens1[3]
		ALT= tokens1[4]
		QUAL= tokens1[5]
		FILTER= tokens1[6]
		INFO= tokens1[7]
		FORMAT= tokens1[8]
		indivs = tokens1[9]

		#1. Skip non-RS
		if not(keepnovel):
			ID = ID.strip()

			if not(ID=='.' or ID[0:2]=="no"):
				continue

			#2. Left with new variants. Now we filter the filter column
			
		try:
		##Check 1a: split the FORMAT column, look for AL identifier position
		#Do this every row? Yes, sadly the format might change between rows.
			AL_index = -1
			formats = FORMAT.split(':')
			for v in xrange(0,len(formats)):
				if formats[v] == AL:
					AL_index = v;
					break;
			
			if(AL_index==-1):
				print >> sys.stderr, "Unable to find AL index in FORMAT column"
				exit(-1);
						
						
		## Check 1b: Check the AL position in the individual
			effect = indivs.split(':')[AL_index]
			#effects = effects.split(',')
			
			# Scan all effects for matching effect. If just one found, keep
			# Each effect index can specify multiple effects
			effect_found = False
			uneffect_found = False

			# Nothing specified
			if filter_name==[] and unfilter_name==[]:
				effect_found = True	
				uneffect_found = False

			else:
				eff_AND = 0			# increments for each matching effect specified with + 
				uneff_AND = 0		# increments for each matching uneffect specified with +
				## Filter scenarios
			#	for effect in effects:

				if filter_name!=[]:

					for fil in filter_name:
						concacs = fil.split('+')

						if len(concacs)>1:
							for c in concacs:
								if re.findall('(('+c+'))', effect):
									eff_AND += 1

							if eff_AND==len(concacs):
								effect_found = True

						else:
							if re.findall('(('+fil+'))', effect):
								effect_found = True

				if unfilter_name!=[]:

					for unfil in unfilter_name:
						concacs = unfil.split('+')

						if len(concacs)>1:
							for c in concacs:
								if re.findall('(('+c+'))', effect):
									uneff_AND += 1
							if uneff_AND==len(concacs):
								uneffect_found = True

						else:
							if re.findall('(('+unfil+'))', effect):
								uneffect_found = True


#				print effects, effect_found, uneffect_found

#				if filter_name==[] and unfilter_name==[]:
#					#Keep line
				if filter_name!=[] and unfilter_name==[]:
					if not(effect_found):
						continue

				elif unfilter_name!=[] and filter_name==[]:
					if not(effect_found) and uneffect_found:
#						print effects, effect_found, uneffect_found
						continue

				elif unfilter_name!=[] and unfilter_name!=[]:
					if effect_found and uneffect_found:
						continue
#					elif effect_found and not(uneffect_found):
#						# Keep line
					elif not(effect_found) and uneffect_found:
						continue
					elif not(effect_found) and not(uneffect_found):
						continue

		except IndexError:
			print line
			exit(-1)

		# If found, either print line, or standardize output
		line = line.splitlines()[0]

		## Now we parse the INFO columns
		tokens = INFO.split(";")

		depth_index = -1
        
		for i in xrange(len(tokens)):
			if tokens[i].startswith("DP="):
				depth_index = i
				break
        
#		print depth_index, tokens[depth_index]

		depth = -1
		try:
			depth = int(tokens[depth_index].split("=")[1])
		except ValueError:
			print >> sys.stderr, "[Error] Could not parse %s" % tokens[depth_index]
		except IndexError:
			print >> sys.stderr, "[Error] Could not find index %s" % line
	
		if depth < depth_lim:
			continue #Skip line

		### PRINT LINE
		if standard:
			variant_lines.append(line+'\t'+str(het))
		else:
			print line

	# OUTPUT accumulated data
	if standard:
		printStandardOut( variant_lines )

	if len(questionable_lines)>1:
		fw = openFile("indels.txt", 'a')
		fw.write("Questionable lines: %s\n" % file)
		fw.write(reduce(lambda x,y: str(x)+'\n'+str(y), questionable_lines))
		fw.write("\n\n")
		fw.close()
		print >> sys.stderr, "\n%d insertions/deletions written to indels.txt" % len(questionable_lines)
		

try:
	files, filt, unfilt, depth, standard, keepnovel = parse_args()

	for file in files:
		process(file, filt, unfilt, depth, standard, keepnovel)

	
except KeyboardInterrupt:
	print >> sys.stderr, "\r[Terminated]"
	

