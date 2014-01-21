#!/usr/bin/env python

import re, sys

name= sys.argv[0]

#0. A R G U M E N T    H A N D L I N G
def usage():
	print >> sys.stdout, '''

      |\ | | | | | | | | | | /|
      || | | | | | | | | | | ||
      ||_|_|_|_|_|_|_|_|_|_|_||
      \____________ver. 2.0.1__/
   .-.. ..-..-..  ..-. .-..-..  ..-. 
   |  |-||( | ||\/||-  |  | ||\/||(  
   `-'' `' '`-''  ``-' `-'`-''  ``-' 

Usage:  %s <VCF-file> <bed-file> [--asterisk]
  
    and where <bed-file> is specified in BED format:
    <chrN>       <start>     <stop>
e.g. chr2        34037655    44354495

    and --asterisk outputs 1 marker (prepended by an asterisk)
    both just before the coding region and just after.
''' % name
	exit()

def grabInputRanges(inputf):
	input_data = []

	for line in open(inputf,'r'):
		vals = line.split()
		if len(vals)>1:
								# Chr (num only), start, stop
			input_data.append( (vals[0].strip(), int(vals[1]), int(vals[2])) )
		else:
			print >> sys.stdout, "Error on '%s'" % (line)

	return input_data;


def message(string, newline=True):
	if (newline):
		print >> sys.stderr, string
	else:
		print >> sys.stderr, string,
	

def processVCF(wb, ranges, file_o):
	
	def grabHeaders():
		line= ""
		for line in wb:
			print line.splitlines()[0]			# Print headers
			if line.find("POS")!=-1:
				return line;

	
	## Extracts only rows pertaining to the current chromosome

	def grabChromosomeWindow(ranges):
		def chrome(vals):
			return int(vals[0].split("hr")[1])

		# Sort ranges in number order
		ranges=sorted(ranges,key=chrome)

		count=0 ## Count num lines -- for fun


		store_rows=[]
		
		finding=False

		for line in wb:
				line = line.splitlines()[0]
		
				if len(ranges)>0:
					for r in ranges:				
						schrom, sstart, sstop = r;
						
						tokes = line.split('\t')
						fchrom = tokes[0].split('hr')[-1].strip(); 
						schrom = schrom.split('hr')[-1].strip();
						fpos = int(tokes[1])
						
						if (fchrom == schrom ):
							#continue
									
							if ((sstop > fpos) and (fpos > sstart)):
								print line
#								finding=True
#							else:
#								if(finding==True):
#									ranges=ranges[1:]
#								finding=False
							
						if (count%300==0):
							print >> sys.stderr, "\r\t\t\t\t\t", schrom, count,
						
							
						#print >> sys.stderr, len(ranges), count
				count +=1
				
		wb.close();
#		if (len(ranges)!=0):
#			print >> sys.stderr, "\n\nError: Unable to use range:", ranges;
#			exit(-1);
			


	headers = grabHeaders()
	grabChromosomeWindow(ranges)
	

def parse_args():
	lenarg = len(sys.argv)

	files=[];
	asterisk=False
	
	index = lenarg-1

	while index > 0:
		arg = sys.argv[index]

		if arg[0]=='-':
			
			if arg=="--asterisk":
				asterisk=True
				sys.argv.remove(arg)

			else:
				usage()

		index -= 1

	sys.argv.remove(name)
	files = sys.argv # remaining args

	if len(files)!=2:
		usage()

	return files[0], files[1], asterisk




vcf_file, input_file, asterisk = parse_args()

message("Open:%s" % vcf_file.split('/')[-1], False)
try:
	wb=open(vcf_file)
except IOError:
	message("[Error] Could not open the file :(")
	exit()
message("X", False)

message("& Parse ", False)
ranges = grabInputRanges(input_file)
message("X", False)

if len(ranges)==0:
	print >> sys.stderr, "[Error] No ranges specified in input file"
	exit()

vcf_file=vcf_file.rsplit('.',1)[0].replace(' ','_').replace('\\','-')

## VCF files are BIG, so chances are processing needs to be done line-by-line
## rather than reading in everything all at once.
processVCF(wb, ranges, vcf_file);
message("\r\t\t\t\t\t\t\t[DONE]");




