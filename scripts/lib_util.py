
'''lib_util.py 

Set of utilities for manipulating biological data : alignments and trees.
'''

import os.path
import os
import string
import random
import shutil

import utilitaires

spsplitchar = utilitaires.spsplitchar


"""
def pruneTrees(intrees, dis_sp_trees, outfile):
	'''Removes from the intree the taxa specified in the list dis_sp_trees.
	
	Returns the name of the resulting outfile.  	
	'''
	
	treefile=open(intrees, 'r')
	trees=treefile.readlines()
	treefile.close()

	tree=trees[0]
	tmp="tmp_tree_%d"%(random.random()*1000000)
	tmp_tree=open(tmp,'w')
	tmp_tree.writelines(tree)
	tmp_tree.close()
		
	for sp in dis_sp_trees:
		cmd='nopartial %s %s %s > /dev/null'%(tmp, tmp, sp)
		print cmd
		os.system(cmd)
		
	os.rename(tmp, outfile)
	
	return outfile
"""	


def extractSpeciesFromTree(tree_file):
	'''  Fonction prenant en argument un arbre au format phylip et retournant la liste des cles correspondantes.
	'''
	list_sp=[]

	tree_file=open(tree_file, 'r')
	tree=tree_file.readlines()
	tree_file.close()

	tree=tree[0][:-2] # Recup du 1er arbre dans le fichier (sans '\n' et ';')


	tab=string.maketrans('','')
	tree=tree.translate(tab, '()') # Suppression des parentheses

	tree=tree.split(',')

	for x in tree:
		list_sp.append(x.split(':')[0])

	return list_sp	

	
def concatAln(aln_file_list, filename, format="interleaved"):
	'''This function concatenates the alignments given in the filename list.
	
	'''
	
	dico_concat={}
	#print filename
	#print aln_file_list
	for aln_file in aln_file_list:
		#print aln_file
		aln=AlnGenerator(aln_file)
		#print aln
		aln.norm_name()
		dico_aln=aln.get_dico_sp_seq()
		#print dico_aln
		sp_list=aln.get_species()
		#print aln_file
		if len(set(sp_list)) > len(sp_list):
			raise InvalidFile, "Several occurences of the same species in file %s:\n%s"%(aln_file, str(sp_list))
		for sp in sp_list:
			dico_concat[sp] = dico_concat.setdefault(sp, "") + dico_aln[sp]
		#	print "%s : %d,"%(sp, len(dico_concat[sp])),
		#print ''	
	
	#print dico_concat
	aln_concat=AlignmentFromDict(dico_concat, filename, format, sp_list)
	
	if format=='fasta':
		aln_concat.write_fasta(filename=filename)
	elif format=='sequential':
		aln_concat.write_sequential(filename=filename)
	else:
		aln_concat.write_interleaved(filename=filename)	

	#print "#### CONCATALN OKI !!!!!!!!!!!!"
	return aln_concat


class InvalidFile(TypeError):
	"Used to indicate that the infile is of an incorrect type"			
		
class InvalidDict(TypeError):
	"Used to indicate that the dictionnary input is of an incorrect type"

	
class Alignment(object):
	"""Abstract class for generic alignments.  
	"""
	def __init__(self, filename, format):
		self.filename=filename
		self.nb_taxa=0
		self._format=format
		self._length=0
		self._dico_sp_seq={}
		self._list_sp=[]
		self._colwidth=60
	
	def __repr__(self):
		return "Redefinir le comportement de print. "
	
	def get_format(self):
		return self._format
	
	def get_nb_taxa(self):
		return self.nb_taxa
	
	def get_dico_sp_seq(self):
		return self._dico_sp_seq
	
	def get_aln_length(self):
		return self._length

	def get_sequence(self, taxa):
		return self._dico_sp_seq[taxa]

	def add_sequence(self, taxa, seq):
		#print self.get_species()
		if self.get_species().count(taxa)==0:
			self._dico_sp_seq[taxa]=seq
			self._list_sp.append(taxa)
			self.nb_taxa+=1
			print "added sequence %s"%taxa
		else:
			print "taxa %s already in alignment"%taxa
	
	def rm_sequence(self, taxa):
		self._dico_sp_seq.pop(taxa)
		self._list_sp.remove(taxa)
		self.nb_taxa-=1		
		print "removed sequence %s"%taxa
	
	def get_column(self,col):
        	"""Returns a string containing a given column with the taxa order specified in self.get_species()
		
		The nb of column starts at 0. """	
		col_str = ''
		assert col >= 0 and col <= self.get_aln_length()		
		for taxa in self._list_sp:
			col_str += self._dico_sp_seq[taxa][col]
		return col_str
	
	def get_species(self):
		#return self._dico_sp_seq.keys()
		return self._list_sp
	
	def get_colwidth(self):
		return self._colwidth
	
	def set_colwidth(self, colwidth):
		if colwidth>0:
			self._colwidth=colwidth
		else:
			print "Colwidth must be >= 1."
			
	def load_aln(self):
		"""Void function. Specified in inheritants. 
		"""
		print "load_aln of Alignment(object)"
		self._dico_sp_seq={}
	
	def norm_name(self, sep=spsplitchar, field=0):
		'''Normalize the name of sequences by splitting them given the separator, and taking the field indicated.		
		'''
		#print "normalizing sequence names"
		dico_sp_seq=self.get_dico_sp_seq()
		new_dico={}
		new_list_sp=[]
		#for sp in dico_sp_seq.keys():
		for sp in self.get_species():
			#print sp
			new_sp=sp.split(sep)[field]
			new_dico[new_sp]=self._dico_sp_seq[sp]
			new_list_sp.append(new_sp)
		#print new_dico
		self._dico_sp_seq=new_dico
		self._list_sp=new_list_sp
		
		
	def write_fasta(self, filename="", colwidth=60):
		"""Write the alignment into a fasta format. 
		
		Default width of columns set at 60 sites.
		"""
		if not filename:
			fastafile=self.filename+".fasta"
		else:
			fastafile=filename
	
		#print "Written in %s"%fastafile
		
		fastafile=open(fastafile, 'w')
		dico=self.get_dico_sp_seq()
		list_sp=self.get_species()	
		for sp in list_sp:
		#for sp in dico.keys():
			seq=dico[sp]
			lines=[]
			line='>%s\n' %sp
			fastafile.write(line)

			i=0
			while i<len(seq):
				line=seq[i:i+colwidth]+'\n'
				i=i+colwidth
				fastafile.write(line)	
	
		fastafile.close()
		
			
	
	def write_interleaved(self, filename="", colwidth=60):
		"""Write the alignment into an interleaved phylip format. 
		
		Default width of columns set at 60 sites.
		"""
		if not filename:
			phylipfile=self.filename+".phy"
		else:
			phylipfile=filename
		
		#print "Written in %s"%phylip_file
		
		phylipfile=open(phylipfile, 'w')
		
		# Mise en forme du fichier : format phylip "interleaved" :
		# 1eres lignes : nom d'espece apparait :
		phylipfile.write("%d\t%d\n"%(self.get_nb_taxa(), self.get_aln_length()))
		i=0
		dico_sp_seq=self.get_dico_sp_seq()		
		nb_sites=self.get_aln_length()
		list_sp=self.get_species()
		
		#for sp in dico_sp_seq.keys():
		for sp in list_sp:
			line="%s      %s \n" %(sp, dico_sp_seq[sp][i:i+colwidth])
			phylipfile.write(line)
		phylipfile.write('\n')
	
		i=i+colwidth
		
		# Autres lignes : pas besoin de mettre noms d'especes
		while i<nb_sites:
			#for sp in dico_sp_seq.keys():
			for sp in list_sp:
				line="%s \n" %(dico_sp_seq[sp][i:i+colwidth])
				phylipfile.write(line)
			phylipfile.write('\n')
			i=i+colwidth	
		
		phylipfile.close()
		
	
	def write_sequential(self, filename="", colwidth=60):
		"""Write the alignment into a sequential phylip format. 
		
		"""
		if not filename:
			phylipfile=self.filename+".phy"
		else:
			phylipfile=filename
			
		#print "Written in %s"%phylipfile
		
		phylipfile=open(phylipfile, 'w')	
			
		dico_sp_seq=self.get_dico_sp_seq()
		list_sp=self.get_species()
		
		# Mise en forme du fichier : format phylip "sequential" :
		phylipfile.write("%d\t%d\n"%(self.get_nb_taxa(), self.get_aln_length()))		
		#for sp in dico_sp_seq.keys():
		for sp in list_sp:
			line="%s      %s \n" %(sp, dico_sp_seq[sp])
			phylipfile.write(line)	
		
		phylipfile.close()	





class AlignmentFromDict(Alignment):
	'''Class of aln that relying on a dico of species name as keys with related sequences as values.
	
	The inherited constructor is redefined.
	'''
	
	def __init__(self, dico, filename, format, sp_list):
		self.filename=filename
		self._format=format
		self._colwidth=60
		self.nb_taxa=0
		self._length=0
		self._dico_sp_seq=self.load_aln(dico)		
		if sp_list:
			self._list_sp=sp_list
		else:
			self._list_sp=self._dico_sp_seq.keys()
	
	def __repr__(self):
		format=self.get_format()
		dico=self.get_dico_sp_seq()		
		colwidth=self.get_colwidth()
		list_sp=self.get_species()
		lines=""
		
		if format=='fasta':
			for sp in list_sp:
			#for sp in dico.keys():
				seq=dico[sp]
				lines+='>%s\n'%sp								

				i=0
				while i<len(seq):
					lines+=seq[i:i+colwidth]
					lines+='\n'
					i=i+colwidth	
			return lines
			
		else:
			nb_sites=self.get_aln_length()

			lines=""
			lines+='%d\t%d\n'%(self.get_nb_taxa(), nb_sites)		
			# 1eres lignes : nom d'espece apparait :
			i=0
			#for sp in dico.keys():
			for sp in list_sp:
				lines+="%s      %s \n" %(sp, dico[sp][i:i+colwidth])
			lines+='\n'
		
			i=i+colwidth		
			# Autres lignes : pas besoin de mettre noms d'especes
			while i<nb_sites:
				for sp in list_sp:
				#for sp in dico.keys():
					lines+="%s \n" %(dico[sp][i:i+colwidth])
				lines+='\n'
				i=i+colwidth	

			return lines
		
	
	def load_aln(self, dico):
		"""Check wether the underlying dictionnary is coherent. 
		
		"""
		#print dico
		length_aln=int(len(dico[dico.keys()[0]]))
		if self._format is not 'fasta':
			for el in dico:
				#print "##### Longueur de la sequence %s : %d VS longueur de l'aln : %d"%(el, len(dico[el]), length_aln)
				if not len(dico[el])==length_aln:
					raise InvalidDict, "The specified dictionnary %s has not the standard format for an alignment.\nSequence %s : %d VS length of 1st sequence: %d"%(self.filename,el, len(dico[el]), length_aln)
				
		self.nb_taxa=len(dico.keys())
		self._length=length_aln
				
		return dico

	

				

class AlnFasta(Alignment):
	"""Alignment class for Fasta format.
	"""

	def __init__(self, filename, format):
		"""Constructor.
		"""
		Alignment.__init__(self, filename, format)
		self.load_aln()
		#print "%s in %s format : %d taxons, %d sites"%(self.filename, self._format, self.nb_taxa, self._length)
	
	def __repr__(self):
		dico=self.get_dico_sp_seq()		
		colwidth=self.get_colwidth()
		list_sp=self.get_species()
		lines=""
		#for sp in dico.keys():
		for sp in list_sp:
			seq=dico[sp]
			lines+='>%s\n'%sp

			i=0
			while i<len(seq):
				lines+=seq[i:i+colwidth]
				lines+='\n'
				i=i+colwidth	
		return lines
		
	def load_aln(self):
		"""Initializing function to load data alignment.
		
		Specific of the fasta format
		"""
		#print "load_aln de AlnFasta(Alignment)"
		infile=open(self.filename,'r')
		inlines=infile.readlines()
		infile.close()
	
		#nb_sp=int(inlines[0].split()[0])
		
		#self._list_sp=[]
		i=0
		for line in inlines:
			if line !='\n':
				if line.startswith(">"):
					#print line
					i+=1
					#taxa=utilitaires.mvWhiteSpaces(line)[1:11]
					# Sans limitation de taille du champ "taxon"
					taxa=line.strip(">\n").rsplit(spsplitchar, 1)[0] # supprime deja le '\n' et le '>' et separe espece et nom de gene !!! string 'spsplitchar' est a changer en fonction de la nature des donnees : CDSs -> "." ; proteines -> "_"
					#~ #print string.join(f)
					#~ taxa=f[0].rsplit(spsplitchar, 1)[0]        # remove trailing replicon identifier
					self._list_sp.append(taxa)
					self._dico_sp_seq[taxa]=""
				else:
					self._dico_sp_seq[taxa]+=utilitaires.mvWhiteSpaces(line[:-1]) 	
				#i+=1
		#print self._dico_sp_seq
		self.nb_taxa=i
		if not self.nb_taxa>0:
			raise InvalidFile, "The specified fasta file %s has no species field. "%(self.filename)
		else:
			
			lens=[]
			for k in self._dico_sp_seq.values():
				lens.append(len(k))
			if len(lens)>0:	
				self._length=max(lens)
			else:
				print "Warning ! The provided fasta alignment %s has 0 sites."%(self.filename)
			#print lens

class AlnPhylip(Alignment):
	"""Alignment class for Phylip format.
	"""
	def __init__(self, filename, format):
		"""Constructor
		"""
		Alignment.__init__(self, filename, format)
		self.load_aln()
		#print "%s in %s format : %d taxons, %d sites"%(self.filename, self._format, self.nb_taxa, self._length)
	
	
	def __repr__(self):
		
		# Mise en forme du fichier : format phylip "interleaved" :
		
		dico_sp_seq=self.get_dico_sp_seq()		
		colwidth=self.get_colwidth()
		nb_sites=self.get_aln_length()
		list_sp=self.get_species()
		
		lines=""
		lines+='%d\t%d\n'%(self.get_nb_taxa(), nb_sites)		
		# 1eres lignes : nom d'espece apparait :
		i=0
		#for sp in dico_sp_seq.keys():
		for sp in list_sp:
			lines+="%s      %s \n" %(sp, dico_sp_seq[sp][i:i+colwidth])
		lines+='\n'
	
		i=i+colwidth		
		# Autres lignes : pas besoin de mettre noms d'especes
		while i<nb_sites:
			#for sp in dico_sp_seq.keys():
			for sp in list_sp:
				lines+="%s \n" %(dico_sp_seq[sp][i:i+colwidth])
			lines+='\n'
			i=i+colwidth	

		return lines
	
	
	
	def load_aln(self):
		"""Initializing function to load data alignment.
		
		Specific of the phylip format
		"""
		#print "load_aln de AlnPhylip(Alignment)"
		
		infile=open(self.filename,'r')
		inlines=infile.readlines()
		infile.close()
	
		self.nb_taxa=int(inlines[0].split()[0])
		self._length=int(inlines[0].split()[1])
		
		if not self.nb_taxa>0:
			raise InvalidFile, "The specified phylip file %s has no species. "%(self.filename)
		
		#self._list_sp=[]
		if self._length>0:
			if self._format == "interleaved":
				i=0
				for line in inlines[1:]:
					if line !='\n':
						fields=line.split()
						if len(fields) == 2:
							sp=fields[0]
							seq=utilitaires.mvWhiteSpaces(fields[1])		
							self._list_sp.append(sp)
							#print sp
							self._dico_sp_seq[sp]=seq
							
						if len(fields) == 1:
							sp=self._list_sp[i%self.nb_taxa]
							#print sp
							seq=utilitaires.mvWhiteSpaces(fields[0])
							self._dico_sp_seq[sp]+=seq
						i+=1
						
			if self._format == "sequential":	
				i=0
				for line in inlines[1:]:
					if line !='\n':
						fields=line.split()
						if len(fields) == 2:
							sp=fields[0]
							seq=utilitaires.mvWhiteSpaces(fields[1])		
							self._list_sp.append(sp)
							self._dico_sp_seq[sp]=seq
							i+=1
				if i!= self.nb_taxa:
					raise InvalidFile, "The specified phylip sequential file %s has not the standard format. "%(self.filename)	
		else:
			print "Warning ! The provided phylip alignment %s has 0 sites."%(self.filename)
		
def AlnGenerator(filename):
	"""Function creating the appropriate Alignment object according to the alignment format. 
	
	Returns an object AlnFasta or AlnPhylip. 
	"""
	fich=open(filename,'r')
	firstline=fich.readline()
	
	format=""
	
	if firstline.startswith(">"):
		format="fasta"
	else:
		fields=firstline.split()
		if len(fields) >= 2:
			if fields[0].isdigit() and fields[1].isdigit():
				format="phylip"
				i=0
				line=fich.readline()
				while line:
				
					if line!='\n': 
						if len(line.split())==1:
							format="interleaved"
							break
						if len(line.split())==2:
							format="sequential"
					line=fich.readline()		
																					
	fich.close()		
	
	if format == "sequential" or format == "interleaved":
		#print "AlnPhylip"
		return AlnPhylip(filename, format)
	if format == "fasta":
		#print "AlnFasta"
		return AlnFasta(filename, format)
	else: 
		raise InvalidFile, "The specified alignment file %s is of unknown format. "%(filename)
		
		
		
		

