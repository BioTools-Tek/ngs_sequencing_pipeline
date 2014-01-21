#!/usr/bin/env python

import sys
from pprint import pprint

aff_arg="--affecteds"
unaff_arg="--unaffecteds"
sib_arg="--siblings"
nonsib_arg="--nonsiblings"
parents_arg="--parents"

def usage():
	print >> sys.stderr, '''
Print the IDs of a pedfile depending on affectation and sibling

	usage: %s <pedfile.pro> [OPTIONS]
	
	OPTIONS:
	%s
	%s
	%s
	%s
	%s''' % (sys.argv[0],aff_arg,unaff_arg,sib_arg,nonsib_arg, parents_arg)
	
	exit(-1)

affecteds=False
unaffecteds=False
siblings=False
nonsiblings=False
parents=False


if len(sys.argv)<2:
	usage()


pedfile = sys.argv[1];

#Opts
for arg in sys.argv[2:]:
	if arg.startswith(aff_arg):
		affecteds=True
	elif arg.startswith(unaff_arg):
		unaffecteds=True
	elif arg.startswith(sib_arg):
		siblings=True
	elif arg.startswith(nonsib_arg):
		nonsiblings=True
	elif arg.startswith(parents_arg):
		parents=True
	else:
		usage()

#Nothing args given, do everything
if (not unaffecteds) and (not affecteds):
	affecteds=True
	unaffecteds=True
if (not siblings) and (not nonsiblings):
	siblings=True
	nonsiblings=True


class Individual:
	def __init__(self, fid, id, father, mother, gender, affected):
		self.fid = fid
		self.id = id;
		self.gender = gender
		self.affected = affected
		self.father = father
		self.mother = mother

		

class Family:
	def __init__(self, famID, famlines):
		self.famID = famID
		self.lines = famlines
		self.individuals = []		
	
	def makeIndivs(self):
		for line in self.lines:
			tokes = line.split();
			tmp = Individual(int(tokes[0]), int(tokes[1]),int(tokes[2]), int(tokes[3]), int(tokes[4]), int(tokes[5]))
			self.individuals.append(tmp);

	def findAllAffecteds(self,affectation=2): ## Only last generation for recessive
		affids=[]
		for indivs in self.individuals:
			if indivs.affected==affectation:
				affids.append(indivs)
		return affids

	def findAllUnaffecteds(self): ## Only last generation for recessive
		return self.findAllAffecteds(1)

	def filterUnaffecteds(self, array):
		return self.filterAffecteds(array,1)

	def filterAffecteds(self,array,num=2):
		return filter(lambda x: x.affected==num, array)

	def findSibs(self,ids,affectation=-1):
		sibs=[]
		non_existent_parents=[]
		all_ids = getIDs(self.individuals)
		if not ids:
			return sibs
		
		for id1 in ids:
			for id2 in self.individuals:

				if id2.mother not in all_ids:
					non_existent_parents.append(id2.mother)
				if id2.father not in all_ids:
					non_existent_parents.append(id2.father)			

				if(id1.id != id2.id):
					if ((id1.mother == id2.mother) or (id1.father==id2.father)):
						#print >> sys.stderr, id1.mother
						if id1 not in sibs:
							sibs.append(id1)
						if id2 not in sibs:
							sibs.append(id2)
		return sibs


	def findNonSibs(self, ids):
		if not ids:
			return []
		
		fams = self.individuals[:]  #clone
		sibs = self.findSibs(ids)
		
		for s in sibs:
			if s in fams:
				fams.remove(s)

		return fams


def getParents(array):
	tmp_moth=-1
	tmp_fath=-1

	for id in array:
		if tmp_moth==-1:
			tmp_moth= id.mother; tmp_fath= id.father
			continue
		
		if not((id.mother == tmp_moth) and (id.father == tmp_fath)):
#			array.remove(id)
			print >> sys.stderr, "[Error]", id.id, "parents:", id.father, id.mother," compared to ", tmp_fath, tmp_moth
#			exit(-1)
		
	return tmp_fath,tmp_moth

def getIDs(array):
	return map(lambda x: x.id, array)

def printOut(famID, arrayed, parents, aff=True, sib=True):
	array = getIDs(arrayed)

	s=" "
	if array:
		if len(array)==1:
			s=str(array[0])
		else:
			s=reduce(lambda x,y: str(x)+"_"+str(y), array)

	beginning = "Fam:"+str(famID)+'\t'+("Aff" if aff else "Unaff")+' '+("Sibs" if sib else "NonSibs");
	end=""
	if not parents:
		end=':'+s
	else:
		dad, mum = getParents(arrayed)
		end='['+s+"] Parents: "+str(dad)+"+"+str(mum)

	if s!=" ":
		print beginning+end

try:
	f = open(pedfile,'r')
except IOError:
	print >> sys.stderr, "could not find pedfile"
	exit(-1);
	
data = f.readlines();
f.close()

tmp_famid=-1
tmp_lines = []

fam_array = []
i = 0;


while (i < len(data)):

	tmp_lines = []
	fid = 0
	line = ""
	
	while (i < len(data)):
		line = data[i]
		tokes = line.split()
		fid = int(tokes[0].strip())
	
		if (tmp_famid==-1):
			tmp_famid = fid

		if fid!=tmp_famid:
			break
		
		tmp_lines.append(data[i])
		i +=1
	
	f = Family(tmp_famid, tmp_lines)
	tmp_famid = fid
	f.makeIndivs()
	
	if affecteds:
		indiv_arr = f.findAllAffecteds()

		if siblings:
			aff_sibs = f.filterAffecteds(f.findSibs(indiv_arr))
			printOut(f.famID, aff_sibs, parents, aff=True, sib=True)
		if nonsiblings:
			aff_nsibs = f.filterAffecteds(f.findNonSibs(indiv_arr))
			printOut(f.famID, aff_nsibs, parents, aff=True, sib=False)
	
	
	if unaffecteds:
		indiv_arr = f.findAllUnaffecteds()
		
		if siblings:
			unaff_sibs= f.filterUnaffecteds(f.findSibs(indiv_arr))
			printOut(f.famID, unaff_sibs, parents, aff=False, sib=True)

		if nonsiblings:
			unaff_nsibs = f.filterUnaffecteds(f.findNonSibs(indiv_arr))
			printOut(f.famID, unaff_nsibs, parents, aff=False, sib=False)
		
	
