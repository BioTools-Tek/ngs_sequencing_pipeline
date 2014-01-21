#!/usr/bin/env python

#mct
import sys

name=sys.argv[0]
delim="~~~~"

def usage():
	print >> sys.stderr, '''
ver 0.1
Kicks out uninformative variants from an affected if its unaffected sibling has the same zygosity

Usage: %s <affected> <unaffected> --zygid=ZYG (recessive|dominant) > out.txt
    
    --zygid=ZYG             Column denoting zygosity
    (recessive|dominant)    kick variant if both are HOM for recessive.
                         OR kick variant if unaff is HET for dominant. (Default)
        
    The file should have no "HETEROZY?" lines

''' % name
	exit()

def parseArgs():
	if len(sys.argv)!=5 :
		usage()

	zygid = sys.argv[3].split("--zygid=")
	zygid = zygid[1] if (not len(zygid)!=2) else usage()
	
	aff_file=sys.argv[1]
	unaff_file=sys.argv[2]
	
	recessive=-1
	
	if (sys.argv[4].strip().startswith("recessive")):
		recessive=True
	elif (sys.argv[4].strip().startswith("dominant")):
		recessive=False
	else:
		print >> sys.stderr, "Error: Need to give 'recessive' or 'dominant'"
		exit(-1)
	
	print >> sys.stderr, '-'*40
	print >> sys.stderr, "Looking for","recessive" if recessive else "dominant", "variants only. ", "Using",zygid,"for zygosity markers\n"

	return aff_file, unaff_file, zygid, recessive


class FileData:
	
	def __init__(self, filename):
		self.filename = filename
		self.headers = []
		self.linemap = {}	# uniq_pos, corresponding_line
		self.zygmap = {}	# uniq_pos, zygmap

	def retrieveLine(self, external_cp):
		try:
			return self.linemap[external_cp]
		except KeyError:
			return -1
			
	def retrieveExists(self,external_cp):
		return (external_cp in self.linemap)
				

def makeChromPosZyg(string, zygid):
	vals = string.split('\t')
	if(len(vals)>1):
		#find zygosity column
		format = vals[-2].split(':')
		try:
			zy_index = format.index(zygid)
		except ValueError:
			print >> sys.stderr, "Cannot find", zygid, "for\n", line
			exit(-1)
			
		values = vals[-1].split(':')
		try:
			zygosity = values[zy_index].split(',')[-1].strip()
		except IndexError:
			print >> sys.stderr, "Could not find zygosity in", values
			exit(-1)

		position=-1
		try:
			position=int(vals[1])
		except ValueError:
			return -1

		if string[0]=='*':
			return -1

		if vals[0].startswith("chr"):
			chrome = vals[0][3:]
			return ((chrome,position), zygosity)

	return -1


def grabDataFrom(filename, zygid):
	count = 0
	f = FileData(filename)
	
	read = openFile(filename)
	print >> sys.stderr, "UnAffected:", filename.split('/')[-1],
	
	for line in read:
		line = line.splitlines()[0]
		
		if (line.startswith('#')):
			f.headers.append(line);
		else:
			count += 1;
			print >> sys.stderr, "\r\t\t\t\t#:%d        " % count,
			cp, zyg = makeChromPosZyg(line, zygid)

			if(cp!=-1):
				f.linemap[cp] = line;
				f.zygmap[cp] = zyg;


	print >> sys.stderr, ""
	return f


def removeDataFrom(filename, zygid, unaff_data, rejects, recessive):
	count = 0;
	rej_count=0;
	common_count = 0;
	notin_count = 0;
	check_data = unaff_data.zygmap.keys()
	
	read = openFile(filename)
	
	print >> sys.stderr, "Affected:", filename.split('/')[-1],
	
	for line in read:
		line = line.splitlines()[0]
		
		if (line.startswith('#')):
			print line
			continue
		
		count += 1
		print >> sys.stderr, "\r\t\t\t\t#:%d        " % count,
		cp, zyg = makeChromPosZyg(line, zygid)
		
		if cp==-1:
			print >> sys.stderr,'\n', line
		else:
			#Check variant exists
			if cp in unaff_data.zygmap:
				unaff_zyg = unaff_data.zygmap[cp]
				
				# File should not contain any heterozy lines, but just in case -- strip off the question marks
				unaff_zyg = unaff_zyg.split('?')[0].strip()
				zyg = zyg.split('?')[0].strip()
				
				# If recessive -- affected is HOM. Kick variant if unaffected is also HOM
				# If dominant -- affected is HET or HOM. Kick variant if unaffected is HET
				
				#Possible scenarios:
				#	Aff 	|	Unaff	|
				#========================
				#	HET 	|	HET 	| -- kick if dominant
				#	HOM 	|	HET 	| -- kick if dominant
				
				#	HOM 	|	HOM 	| -- kick if recessive
				#
				# i.e.	kick if equal
				#		kick if aff=HOM and unaff=HET
				#
				if recessive:
					if (zyg == unaff_zyg) and (zyg[0:3]=="HOM"):
						print >> rejects, zyg,delim,unaff_zyg,delim,line,delim,unaff_data.linemap[cp]
						rej_count += 1;
						continue
				
				else:
					if unaff_zyg[0:3]=="HET":
						print >> rejects, zyg,delim,unaff_zyg,delim,line,delim,unaff_data.linemap[cp]
						rej_count += 1;
						continue
				
				common_count += 1
				print line
			
			# If variant not in unaffected, print it anyway
			else:
				notin_count += 1;
				print line

				
	print >> sys.stderr, "\nAffected: Common[",common_count,"], NotCommon[",notin_count, "], Rejects[", rej_count,"]"
	print >> sys.stderr, '-'*40

def openFile(file_i, mode='r'):
	try:
		f=open(file_i,mode)
	except IOError:
		print >> sys.stderr, "[Error] Cannot open %s!" % file_i
		exit()
	return f


def main(args):
	affected, unaffected, zygid, recessive = args
	aff = affected.split('/')[-1]; unaff = unaffected.split('/')[-1];
	
	unaffected_data = grabDataFrom(unaffected, zygid)

	rejects = openFile(aff.split('.')[0]+"_rejects.vcf",'w')
	
	#Rejects headers
	print >> rejects, aff,"zygosity",delim,unaff,"zygosity",delim,aff,delim,unaff
	
	removeDataFrom(affected, zygid, unaffected_data, rejects, recessive)


try:
	main(parseArgs());
except KeyboardInterrupt:
	print >> sys.stderr, "\r[Terminated]                  "
	exit()

