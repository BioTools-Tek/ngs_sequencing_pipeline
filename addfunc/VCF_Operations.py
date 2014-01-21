#!/usr/bin/env python

# Common VCF file handling tasks

import sys

version=0.2

"Performs operations across entire VCF files"
class VCFOps:
	ErrorHead = "[VCF Ops Error]: "
	
	#Seperators
	element_sep = ','
	id_sep=':'
	field_sep='\t'

	def __init__(self, filename):
		self.filename = filename
		self.numlines = self.countlines()
		self.field_map = {}
		self.FORMAT_INDEX = -1
		self.IFORMAT_INDEX = -1
		self.FORMAT_IDS_ORDER = []  # map isn't enough, we need order ids appear in VCF
		self.data_linecount = 0;

	def countlines(self):
		i=0;
		with open(self.filename,'r') as f:
			for i, l in enumerate(f):
				pass
		f.close()
		return i + 1


	"Same as above but static for other files"
	@staticmethod
	def countFileLines(filename):
		i=0;
		with open(filename,'r') as f:
			for i, l in enumerate(f):
				pass
		f.close()
		return i + 1


	"Find Id's that have the same length as each other"
	def getSameLengthIds(self):

		if self.IFORMAT_INDEX == -1:
			print >> sys.stderr, VCFOps.ErrorHead, "Please run getFieldIndexes first!"
			exit(-1)
		
		f = open(self.filename,'r')
		for line in f:
			if line[0]=='#':
				continue
			tokens = line.split(VCFOps.field_sep)
			
			data_fields = tokens[self.IFORMAT_INDEX].split(self.id_sep)
			id_fields   = tokens[self.FORMAT_INDEX] .split(self.id_sep)
			
			# Find a line with at least 2 elements in one of the fields, then we can start checking.
			at_least_len2_found = False;
			id_map={}
			
			for index in xrange(len(data_fields)):
				field = data_fields[index].split(self.element_sep)
				length = len( field )
				
				if length not in id_map:
					id_map[length] = []
				
				id_add = id_fields[index].strip();
				if id_add not in ['GT','AD','DP','GQ','PL']: #ignore these fields
					id_map[length].append( id_add )
				
				#If we have two or more fields that have N elements, where N>=2, then this is a good line to examine
				if length>=2 and len(id_map[length])>=2:
					at_least_len2_found = True


			if at_least_len2_found:
				# find key with most number of elements
				key=-1; max=-1
				for k,v in id_map.iteritems():
					len_items = len(v)
					if max < len_items:
						max = len_items
						key = k;
				
				return id_map[key]


	"Updates IFORMAT (data) and FORMAT (ids) indexes, as well as "
	"a map of where each id lies in each of the formats"

	"FORMAT 		ZYG:CH:PCY:AL"
	"IFORMAT		HETERO:c56A>C:pR67Q:Gene1"
	"field_map		{ZYG:0, PCY:2, etc}"
	def getFieldIndexes(self, arg_ids=[]):

		f = open(self.filename,'r')
		for line in f:
			
			if line[0]=='#':
				if line.startswith("#CHROM"):
					self.FORMAT_INDEX = map(lambda x: x.upper(), line.split('\t')).index("FORMAT")
					self.IFORMAT_INDEX = self.FORMAT_INDEX + 1
				continue
			
			tokens = line.split(VCFOps.field_sep)
			
			FORMAT_DATA = tokens[self.FORMAT_INDEX]
			ids = FORMAT_DATA.split(VCFOps.id_sep)
			self.FORMAT_IDS_ORDER = ids;  # update for other functions
			
			if arg_ids[0]!='':              # check argument ids for existence (if given)
				for a_id in arg_ids:
					if a_id not in ids:
						print >> sys.stderr, VCFOps.ErrorHead, "could not find field '%s' in data!" % a_id
						exit(-1)
			
			for id in ids:                # map all ids
				self.field_map[id] = ids.index(id)
			
			break # everything found, break and close file
		f.close()

		if ( -1 in (self.IFORMAT_INDEX, self.FORMAT_INDEX)):
			print >> sys.stderr, VCFOps.ErrorHead, "Could not find IFORMAT or FORMAT indexes"
			exit(-1)


	"Asserts that the number of elements in each field given by the map "
	"have the same number of elements delimited by splitter, and returns"
	"the data as a map using the same keys as the argument map"
	def getElementsFromData(self, array_of_data, check):
		self.data_linecount +=1
		
		lengths=-1
		res_map={}
		
		
		if check[0]=='':
			print >> sys.stderr, "check list is empty!"
			exit(-1)
		
		for name in check:
			index = self.field_map[name]
			
			data = map(lambda x: x.strip(), array_of_data[index].split(VCFOps.element_sep))
			len_data = len(data)
			if  lengths == -1:
				lengths = len_data
			else:
				if len_data != lengths:
					print >> sys.stderr, VCFOps.ErrorHead, "Line:", self.data_linecount, " - Length of", name, "does not match other fields!"
					print >> sys.stderr, array_of_data, "   ==", name, index, len_data, lengths, check
					exit(-1)
			res_map[name]=data
					
		return res_map
	
	
	
