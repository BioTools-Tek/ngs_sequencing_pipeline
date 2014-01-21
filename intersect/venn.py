#!/usr/bin/env python

#mct

## Searches for common positions between different input files of the same
## coding window. It produces a venn diagram of the intermediary positions
##
import sys

name=sys.argv[0]

def usage():
	print >> sys.stderr, '''
ver 2.3
Compares the overlap of positions between two or more files

Usage: %s <[VCF FILES]> [OPTIONS] > out.txt
              
       --chr NN     compare only one particular chromosome (e.g 01)
       --columns    print the lines from each file too at each position
       --overlap    print only the overlapping lines
       --sep        seperator used for printing out lines across
                    different files. Default is "::"
       --printTot   print totals. Default is false.

''' % name
	exit()


def parseArgs():
     lenarg = len(sys.argv)

     chrome=-1  #defaults
     columns=False
     overlap=False
     printTot=False
     sep="::"
     files=[]

     index=lenarg-1
     while index > 0:
        arg = sys.argv[index]

        if arg.startswith('-'):

           if arg=="--columns":
               columns=True
               sys.argv.remove(arg)
               
           elif arg=="--overlap":
               overlap=True
               sys.argv.remove(arg)

           elif arg=="--printTot":
               printTot=True
               sys.argv.remove(arg)

           elif arg=="--chr":
               try:
                  chrome=sys.argv[index+1]
                  sys.argv.remove(arg)
                  sys.argv.remove(chrome)
                  chrome=int(chrome)
                  if chrome > 22 or chrome==0:
                      print >> sys.stderr, "[Error] Invalid Chromosome"
                      exit(-1)
               except ValueError:
                  print >> sys.stderr, "[Error] Not a Chromosome"
                  exit(-1)
               except IndexError:
                  print >> sys.stderr, "[Error] Must specify a chromosome!"     
                  exit(-1)
                  
           elif arg=="--sep":
               try:
                  sep=sys.argv[index+1].strip()
                  sys.argv.remove(arg)
                  sys.argv.remove(sep)
                  if(len(sep)==0):
					  print >> sys.stderr, "[Error] Bad seperator"
					  exit(-1)                  
               except ValueError:
                  print >> sys.stderr, "[Error] Bad seperator"
                  exit(-1)
               except IndexError:
                  print >> sys.stderr, "[Error] Must specify a seperator!"     
                  exit(-1)


           else:
               usage()

        index -= 1

     
     sys.argv.remove(name)
     files=sys.argv

     if len(files)<2:
         usage()

     return files, chrome, columns, overlap, sep, printTot


class FileData:
	
	def __init__(self, filename):
		self.filename = filename
		self.headers = []
		self.linemap = {}

	def retrieveLine(self, external_cp):
		try:
			return self.linemap[external_cp]
		except KeyError:
			return -1
			
	def retrieveExists(self,external_cp):
		return (external_cp in self.linemap)
				

def makeChromPos(string):
	vals = string.split()
	if(len(vals)>1):
		position=-1
		try:
			position=int(vals[1])
		except ValueError:
			return -1

		if string[0]=='*':
			return -1

		if vals[0].startswith("chr"):
			chrome = vals[0][3:]
			return (chrome,position)
		else:
			print >> sys.stderr, "FLAG:", chrome, position
			return -1

	return -1



def populateUnique_storeFile(read, filename, populate_array, chrome=-1):

	count = 0
	f = FileData(filename)

	for line in read:
		line = line.splitlines()[0]
		
		if (line.startswith('#')):
			f.headers.append(line);
		else:
			print >> sys.stderr, "\r\t\t\t\t\t\t#:%d        " % count,
			uniq = makeChromPos(line)

			if(uniq!=-1):
				f.linemap[uniq] = line;
				
				chro, poso = uniq

				if chro==chrome or chrome==-1 :
					if not(uniq in populate_array):
						populate_array.append(uniq)

		count += 1;

	print >> sys.stderr, ""
	return f


def openFile(file_i, write=False):
	try:
		if write:
			f=open(file_i,'w')
		else:
			f=open(file_i,'r')
	except IOError:
		print >> sys.stderr, "[Error] Cannot open %s!" % file_i
		exit()
	return f
	

#Check for duplicates as you print out
def printCommon(fd_list, uniq, chromosome, columns, overlap, sep, pTot):	

	def printExistingHeaders():
		maxlen=-1
		header_array=[]
		
		#Find max len
		for ff in fd_list:
			newlen = len(ff.headers)
			maxlen = newlen if (maxlen<newlen) else maxlen
			
		#Populate headers column by column
		for index in xrange(-1,maxlen):			
			tmp_rows=[]
			for file in fd_list:
				#Create a dummy first column
#				tmp_rows.append(sep)

				try:
					if index==-1:
						head = file.filename.strip();
					else:
						head = file.headers[index].strip()
				except IndexError:
					head = " "
				tmp_rows.append(head)
			
			header_array.append(tmp_rows);

		#Print
		for h in header_array:
			sys.stdout.write(sep)
			for i in h:
				sys.stdout.write("%s%s%s" % (sep,i,sep))
			sys.stdout.write("%s\n" % sep)


	def makeHeaders():
		res="Chr\tPhysPos\t\t"
		for title in fd_list:
			res += title.filename+"\t"

		res += "\tSum\t\t"

		if columns:
			res += sep
			for title in fd_list:
				res += sep+title.filename+"_ROWS"+sep
			res += sep
		
		if columns:
			printExistingHeaders()

		print res

				

	makeHeaders();
		
	sum_columns={}			## Number of markers appearing in both file and uniq
	uniq_columns={}			## Number of unique markers appearing in one file only

	uniq_matches = 0
	for files in fd_list:
		sum_columns[files.filename] = 0
		uniq_columns[files.filename] = 0

	## Scan rows
	for ext_tuple in uniq:
		u_chrome, u_pos = ext_tuple		##  Process all chromosomes
			
		if chromosome==-1 or chromosome==u_chrome:
			out_line = "%s\t%d\t\t" % (u_chrome, u_pos)

			sum_row=0
			last_contribute = ""		# Name of last file to add to the current sum
			
			for file in fd_list:
				match = file.retrieveExists(ext_tuple)

				if match:
					sum_columns[file.filename] += 1
					sum_row += 1
					last_contribute = file.filename
					
				match = " " if (match==False) else "1"
				
				out_line += match+"\t"

			## Unique to a single file only	
			if sum_row==1:
				uniq_matches +=1
				uniq_columns[last_contribute] += 1

			## Print out each columns
			data_rows=""
			genes = []
			if columns:
				data_rows=sep
				for file in fd_list:
					row = file.retrieveLine(ext_tuple)
					if row==-1:
						row="-"

					data_rows += sep+row+sep
				data_rows += sep

			out_line += ("\t%d\t\t%s\t" % (sum_row, data_rows))
			
			if overlap:
				if sum_row==len(fd_list):
					print out_line
			else:
				print out_line
			
			


	print  >> sys.stderr, "%s%d%s" % (fillerNChars("",uniq_matches,50," "), uniq_matches, " unique" )
	final_line= ""

	if len(uniq)==0:
		uniq.append(1)

	if pTot:

		# Find largest common set
		small = sys.maxint;
		max_common=""

		for file in fd_list:
			col_vals = sum_columns[file.filename]
			uniq_vals = uniq_columns[file.filename]

			if col_vals < small:
				small = col_vals
				max_common = file.filename

			namer=file.filename.split('/')[-1]


			final_line += (
				"%s:%s%d out of %d\n" % (
						namer,
						fillerNChars(namer+":", str(uniq_vals),50," "),
						uniq_vals,
						col_vals
	#					int( 100 * ( float(col_vals) / float(len(uniq)) ))))
						))


		print >> sys.stderr, final_line
		print >> sys.stderr, "Common set:%s%s" % (
				fillerNChars("Common set:", max_common, 56, " "), max_common)


def populate_comparison_arrays(file_list, chrome=-1):
	unique_pos=[]
	fd_array = []

	## Not complete, perform a quick input file wide search for common positions in all
	for file_i in file_list:
		try:
			print >> sys.stderr, "\r Started: %s   " % file_i.split('/')[-1], 
			f = open(file_i,'r')
		except IOError:
			print >>  sys.stderr, "[Error] File %s not found" % file_i
			inp = raw_input("Ignore?[y/n]")

			if inp[0]=="n" or inp[0]=="N":
				print >>  sys.stderr, "\r[Terminated]             "
				exit()
			else:
				continue
			
		file_data = populateUnique_storeFile(f, file_i, unique_pos, chrome)
		f.close()

		fd_array.append(file_data)

	unique_pos.sort(key=lambda x: ( (int(x[0]) if x[0].isdigit() else x[0]),int(x[1])))

	return fd_array, unique_pos


def fillerNChars(front, back, n, char):
	num_chars1=len(str(front))
	num_chars2=len(str(back))
	
	filler = n - (num_chars1 + num_chars2)

	return char*filler


def main(argu):
	files, chromosome, columns, overlap, seperator, pTot = argu

	fd, uniq = populate_comparison_arrays(files,chromosome);
	printCommon(fd,uniq, chromosome, columns, overlap, seperator, pTot)


try:
	main(parseArgs());
except KeyboardInterrupt:
	print >> sys.stderr, "\r[Terminated]                  "
	exit()


