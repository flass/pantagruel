#! /usr/local/bin/python
# utilitaires.py

import os.path
import os
import string
import random
import shutil
import fpformat
import exceptions
import re

#from Bio import Fasta
# A remettre
#~ spsplitchar = '@'
spsplitchar = '_'
#~ import acnuc	

class EmptyFileException(IOError):
	"Used to indicate that the specified file is empty."


def fileToLines(filename):
	fich=open(filename,'r')
	lines=fich.readlines()
	fich.close()
	
	if len(lines)==0:
		raise EmptyFileException("fileTolLines : file %s empty."%filename)
		#break
		#exit(-1)
	
	'''
	if len(lines)==0:
		print "fileTolLines : file %s empty."%filename	
		exit(0)
	'''	
	return lines

def isBootstrapped(tree):
	'''
	Fonction verifiant si l'arbre en entree (string) est bootstrappe ou non. Renvoie 1 si c'est le cas, 0 sinon
	'''
	dico={'(':'', ',':' ', ':':' ', ')':') '}
	new_tree=translate(dico, tree)
	fields=new_tree.split()

	i=0
	while i <=len(fields)-2:
		f=fields[i]
		if f.count('.')>0:
			f=f.replace('.','')
		if f.isdigit():
			g=fields[i+1]
			if g.count('.')>0:
				g=g.replace('.','')
			if g.isdigit():
				# 2 digits se suivent, on a des bootstraps
				return True
			i+=2	
		else:
			i+=1
	return False


def readAUTreefinderResults(logfile, num_ref_tree, threshold=0.05):
	'''This function read the AU result computed by Treefinder in its log file. 
	
	The threshold is 0.05 by default, and the function returns 1 if the test is accepted, 0 if rejected. 
	'''

	resfile=open(logfile, 'r')
	reslines=resfile.readlines()
	resfile.close()
	nb_tree=0
	
	for line in reslines:                
		if line.startswith(" {PairedSitesSupport"):
			#res_line+="\nTree %d\t"%nb_tree
			res_tests=line.split('->{')[1]
			res_tests=res_tests.split('},')[0]
			res_tests=res_tests.split(',')
			if nb_tree==num_ref_tree:
				au_res=res_tests[5]
        			break
	        	nb_tree+=1
			
	if float(au_res) >= threshold:
		#print "%.5f >= %.5f"%(float(au_res),float(threshold))
		return (1,au_res)
	else:
		#print "%.5f < %.5f"%(float(au_res),float(threshold))
		return (0,au_res)



def rmbp(tree):
	'''
	Fonction suppprimant les valeurs bp dans un arbre string au format newick (finissant avec ';') 
	'''

	tree_without_bp=""
	
	# Attention !!! D'abord verifier que l'arbre comporte bien des bootstraps... 	
	if isBootstrapped(tree):
		f1=tree.split(')')
		for f in f1:		
			s=f.split(':')
			if s[0].replace('.','').isdigit():
				tmp=string.join(s[1:])
				#print "if : '%s'"%tmp
				if not tmp.startswith(' '):
					tmp=tmp.replace(' ',':')
				if tmp.startswith("("):
					tree_without_bp+=tmp				
				elif tmp.startswith(';'):
					tree_without_bp+=")%s"%tmp
				else:	
					tree_without_bp+="):%s"%tmp
				
			else:	
				tmp=string.join(s)
				tmp=tmp.lstrip() # !!! en test
				tmp=tmp.replace(' ',':') # !!! en test
				
				#print "else : '%s'"%tmp
				if tmp.startswith("("):
					tree_without_bp+=tmp				
				elif tmp.startswith(';'):
					tree_without_bp+=")%s"%tmp
				else:
					tree_without_bp+="):%s"%tmp
					
		tree_without_bp=tree_without_bp.replace(' ','')
			
		return tree_without_bp
	else:	
		#print "Pas de modif"
		return tree

def extractbp(tree):
	'''
	Fonction extrayant les valeurs de bp dans un arbre string au format newick (finissant avec ';') 
	'''

	bp=[]
	
	# Attention !!! D'abord verifier que l'arbre comporte bien des bootstraps... 	
	if isBootstrapped(tree):
		#regex=re.compile('\)[0-9]+\.?[0-9]*:')
		#regex=re.compile('\)[0-9]+\.(?<=[0-9]*):')
		regex=re.compile('\)[0-9]+\.?[0-9]*:')
		res=regex.findall(tree)	
		#res=regex.finditer(tree)	
		for r in res:
			bp.append(r[1:-1])
			
		return bp
	else:	
		#print "Pas de modif"
		return bp





def scientificNotationToFloat(s):
	'''
	Utilitaire prenant en argument une chaine de caractere representant une valeur numerique en ecriture scientifique. 
	Retourne cette valeur sous forme d'une chaine de caractere d'un flottant (ecriture decimale). 
	'''
	f1=s.split('.')
	un=f1[0]
	tmp=f1[1]
	f2=tmp.split('e')
	deux=f2[0]
	trois=f2[1]	
	nb_dig=len(deux)+abs(int(trois))
	f=fpformat.fix(s, nb_dig)
	return f




def get_randint():
	rand=random.random()*1000000
	rand_str=str(rand).split('.')[0]
	if len(rand_str)!=6:
		diff=(len(rand_str)-6)
		if diff<0:
			#print "diff<0 : %d"%diff
			i=0
			zeros=""
			while i < (6-len(rand_str)):
				zeros+="0"
				i+=1
			randname="%s%s"%(zeros,rand_str)
		else:
			randname=rand_str[0:6]
	else:
		randname=rand_str
	#print rand
	#print randname		
	return int(randname)


def getGBFilenames(wd_fam, phylum_norm, fam_list):
	gb_list=[]
	for fam in fam_list:	
		gb_file="%s_%s.fasta-aln-gb"%(phylum_norm.upper(),fam)
		gb_file="%s/%s"%(wd_fam,gb_file)
		if os.path.exists(gb_file):
			gb_list.append(gb_file)

	return gb_list

def getFastaAlnFilenames(wd_fam, phylum_norm, fam_list):
	fa_list=[]
	for fam in fam_list:	
		fa_file="%s_%s.fasta-aln"%(phylum_norm.upper(),fam)
		#fa_file="%s_%s.fasta-aln"%(phylum_norm,fam)	
		fa_file="%s/%s"%(wd_fam,fa_file)
		if os.path.exists(fa_file):
			fa_list.append(fa_file)

	return fa_list

def getPhylipAlnFilenames(wd_fam, phylum_norm, fam_list):
	fa_list=[]
	for fam in fam_list:	
		fa_file="%s_%s.fasta-aln.phy"%(phylum_norm.upper(),fam)
		#fa_file="%s_%s.fasta-aln"%(phylum_norm,fam)	
		fa_file="%s/%s"%(wd_fam,fa_file)
		if os.path.exists(fa_file):
			fa_list.append(fa_file)

	return fa_list

def getPhylipAlnGBFilenames(wd_fam, phylum_norm, fam_list):
	fa_list=[]
	for fam in fam_list:	
		fa_file="%s_%s.fasta-aln-gb.phy"%(phylum_norm.upper(),fam)
		#fa_file="%s_%s.fasta-aln"%(phylum_norm,fam)	
		fa_file="%s/%s"%(wd_fam,fa_file)
		if os.path.exists(fa_file):
			fa_list.append(fa_file)

	return fa_list

def translate(dico_eq, string):
	'''
	Translate the words given as key of the input dictionnary into the word of the corresponding field.
	'''
		
	new_str=string
	for key in dico_eq:
		#print key
		#print dico_eq[key]
		new_str=new_str.replace(key,dico_eq[key])
	#print new_str	
	return new_str


def fastaConcat(fasta_infile, concat_file):
	'''
	 Fonction faisant la concatenation des fichiers fasta donnes en liste dans le 1er argument. 
	 Le 2eme argument est le fichier concatenat sortie. 	 
	'''

	if os.path.exists(fasta_infile):
		infile=open(fasta_infile,'r')
		fasta_list_lines=infile.readlines()
		infile.close()
		
		fasta_list_files=[]
		# Creation de la liste de fichiers fasta
		for line in fasta_list_lines:
			fasta_list_files.append(line[:-1])
		
		# Initialisation du concatenat : 
		if len(fasta_list_files>=2):	
			file1=fasta_list_files[0]
			file2=fasta_list_files[1]
	
			fasta_concat2files(file1, file2, concat_file, colwidth)
			i=2
			while i < len(fasta_file_list):
				file2=fasta_file_list[i]
				fasta_concat2files(concat_file, file2, concat_file, colwidth)
				i+=1	
		else:
		
			print "The list of fasta files is badly initialized. No more than 1 file ? "
			return 0
	else:
		print "The specified file doesn't exist."	
		return 0





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
		#print sp
		cmd='nopartial %s %s %s > /dev/null'%(tmp, tmp, sp)
		#print cmd
		os.system(cmd)
		
	os.rename(tmp, outfile)
	
	return outfile
	

def getMinRankMin(liste):
	''' Fonction retournant la valeur minimale et l'indice du min d'une liste de valeurs numeriques.	
	
	Rem : existe deja dans python (fonction min())
	'''
	rank=0
	mini=liste[0]
	for i in range(len(liste)):
		el=liste[i]
		if el<mini:
			mini=el
			rank=i
	return (mini,rank)		
	
	
def FastaToPhylip(fasta_file):
	'''Fonction de conversion de fichiers Fasta au format Phylip (interleaved)

	Renvoie le nom du fichier Phylip
	'''

	aln=lib_util.AlnGenerator(fasta)
	aln.write_interleaved()
	return "%s.phy"%fasta	



#*******************************************
# Fonction de conversion de fichiers Fasta au format Phylip (interleaved)
#*******************************************
def fastaToPhylipOLD(fasta_file, colwidth=50):
	'''Fonction de conversion de fichiers Fasta au format Phylip (interleaved)
	'''

	phylip_filename=fasta_file+'.phylip'	
	phylip_file=open(phylip_filename, 'w')
	
	fasta_file=open(fasta_file, 'r')
	parser=Fasta.RecordParser()
	iterator=Fasta.Iterator(fasta_file, parser)
	dico_rec={}
	
	
	while 1:
		record=iterator.next()	
		if record is None:
			break
		sp=record.title	
		
		# Normalisation du nom des sequences : nom d'espece en 5 lettres max
		fields_sp=sp.split()
		#print len(fields_sp)
		if len(fields_sp) > 1:
			sp=fields_sp[0]
			#print sp
		fields_sp=sp.split(spsplitchar)		
		if len(fields_sp) > 1:
			#sp=fields_sp[1] # Pour HOGENOM 3
			sp=fields_sp[0] # Pour HOGENOM 4+
			#print sp
		
		while len(sp) < 10:
			sp+=" "
		
		# Suppression des eventuels espaces (si sortie de gblocks par exemple)
		seq=mvWhiteSpaces(record.sequence)
		nb_sites=len(seq) # Nb de sites pour en-tete du fichier phylip
		dico_rec[sp]=seq
		#print nb_sites
	
	nb_taxa=len(dico_rec) # Nb de sequences (taxa) representes dans le fichier pour en-tete du fichier phylip
	print "%d taxa et %d sites" %(nb_taxa, nb_sites)

	# En-tete du fichier :
	phylip_file.write('%d %d \n' %(nb_taxa, nb_sites)) 

	# Mise en forme du fichier : format phylip "interleaved" :
	# 1eres lignes : nom d'espece apparait :
	i=0
	for sp in dico_rec.keys():
		#line=sp+' '+dico_rec[sp][i:i+colwidth]+'\n'
		line="%s %s \n" %(sp, dico_rec[sp][i:i+colwidth])
		phylip_file.write(line)
	phylip_file.write('\n')
	
	i=i+colwidth
	
	# Autres lignes : pas besoin de mettre noms d'especes
	while i<nb_sites:
		for sp in dico_rec.keys():
			#line=dico_rec[sp][i:i+colwidth]+'\n'
			line="%s \n" %(dico_rec[sp][i:i+colwidth])
			phylip_file.write(line)
		phylip_file.write('\n')
		i=i+colwidth
				
	
	fasta_file.close()
	phylip_file.close()

	return phylip_filename
	
#*******************************************
# Fonction qui concatene 2 fichiers au format fasta comportant le meme echantillon taxonomique (nom des sequences du type : QU41332_AGTU0 )
#*******************************************	
def fasta_concat2files(file1, file2, concat_file, colwidth):
	'''Fonction qui concatene 2 fichiers au format fasta comportant le meme echantillon taxonomique (nom des sequences du type : QU41332_AGTU0 )
	'''
	
	#print '* Concatenation de '+file1+' et '+file2

	file1=open(file1, 'r')
	file2=open(file2, 'r')
	colwidth=60
	
	dico1={}
	dico2={}
	dico_res={}
	
	# Extraction des sequences de la 1ere famille au format fasta : 
	parser=Fasta.RecordParser()
	iterator=Fasta.Iterator(file1, parser)
	while 1:	
		record=iterator.next()
		if record is None:
			break
		sp=record.title.split()[0]
		#~ if len(sp.split(spsplitchar)) > 1:
			#~ sp=sp.split(spsplitchar)[1]
		dico1[sp]=mvWhiteSpaces(record.sequence) # Supprime les espaces ajoutes par Gblocks
	
	
	# Extraction des sequences de la 2eme famille au format fasta : 
	parser=Fasta.RecordParser()
	iterator=Fasta.Iterator(file2, parser)	
	while 1:	
		record=iterator.next()
		if record is None:
			break
		sp=record.title.split()[0]
		#~ if len(sp.split(spsplitchar)) > 1:
			#~ sp=sp.split(spsplitchar)[1]
		dico2[sp]=mvWhiteSpaces(record.sequence) # Supprime les espaces ajoutes par Gblocks


	file1.close()	
	file2.close()	
	
	# Concatenation : 
	for sp in dico1.keys():
		seq_cat=dico1[sp]+dico2[sp]
		dico_res[sp]=seq_cat

	
	# Ecriture dans un nouveau fichier :	
	concat_file=open(concat_file, 'w')		
	for sp in dico_res.keys():
		seq=dico_res[sp]
		lines=[]
		line='>%s\n' %sp
		concat_file.write(line)

		i=0
		while i<len(seq):
			line=seq[i:i+colwidth]+'\n'
			i=i+colwidth
			concat_file.write(line)
	
	
	concat_file.close()
	
		
#*******************************************
# Fonction qui :  	
# Lance la concatenation sur une liste de fichiers fasta comportant tous une sequence par especes
# sortie : 1 fichier fasta de nom jack_1_ALPAPROTEO.fasta-aln-gb avec comme identifiant de sequence > AGTU0
# retourne le nom du fichier concatene
#*******************************************	
def fasta_concat(wd, fasta_file_list, phylum, round_name, colwidth=60):
	'''Fonction qui lance la concatenation d'une liste de fichiers fasta comportant tous une sequence par especes

	sortie : 1 fichier fasta de nom wd/round_name.fasta avec comme identifiant de sequence ex: "> AGTU0"
	retourne le nom du fichier concatene
	'''


	#def fasta_concat(wd, fasta_file_list,colwidth=60):
	#print concat_file
	concat_file=getPhylumFilename(wd, round_name, '.fasta')
	
	#print "@@@"
	#print round_name
	#print "@@@"
	
	file1=fasta_file_list[0]
	#print 'file 1 %s' % file1
	file2=fasta_file_list[1]
	#print 'file 2 %s' % file2
	
	fasta_concat2files(file1, file2, concat_file, colwidth)
	i=2
	while i < len(fasta_file_list):
		file2=fasta_file_list[i]
		fasta_concat2files(concat_file, file2, concat_file, colwidth)
		#print 'concat_file %s' % concat_file
		#print 'file 2 %s' % file2
		i+=1
	 
	#print "@@@@@@@@@\n"
	#print fasta_file_list
	
	return concat_file	
	
#*******************************************
# Fonction efectuant un tirage aleatoire sans remise de nb_fam elements de la liste fam_list
# Retourne la liste de ces elements
#*******************************************	
		
def fam_jackknife(fam_list, nb_fam):
	''' Fonction efectuant un tirage aleatoire sans remise de nb_fam elements d'une liste fam_list
	
	Retourne la liste de ces elements
	'''

	# Echantillon de familles prises au hasard (tirage sans remise)
	fam_jack=random.sample(fam_list, nb_fam)
	return fam_jack	
	

	
	

#*******************************************
# Fonction retournant la chaine de caractere en entree sans les espaces
#*******************************************
def mvWhiteSpaces(seq):
	'''Fonction retournant la chaine de caractere en entree sans les espaces
	'''
	
	tab=string.maketrans('','') # table "identite" : cas ou besoin simplement de specifier un caractere a supprimer
	return seq.translate(tab, ' ') # supprime les espaces


#*******************************************
# Fonction retournant un nom de fichier absolu pour un nom de phylum donne de type 'ACTINOBACTERIA (CLASS)'
# entree : repertoire de travail, nom complet du phylum (peut etre compose de 2 mots) et suffixe a ajouter ('.sp')
# sortie : retourne la nom de fichier
#*******************************************
def getPhylumFilename(wd, phylum, suffix):
	'''Fonction retournant un nom de fichier absolu pour un nom de phylum donne de type 'ACTINOBACTERIA (CLASS)'

	entree : repertoire de travail, nom complet du phylum (peut etre compose de 2 mots) et suffixe a ajouter ('.sp')
	sortie : retourne la nom de fichier
	'''

	tab=string.maketrans('/','_')
	phylum=phylum.translate(tab)
	filename=phylum.split()[0]+suffix	
	
	abs_filename=os.path.join(wd, filename)
	#print "'"+wd+"'"
	#print "$"+filename+"$"
	#print abs_filename
	
	return abs_filename



#*******************************************
# Fonction retournant le nom de l'espece (Genre)
#*******************************************
def getGenus(sp_name):
	'''Fonction retournant le nom de l'espece (Genre)
	'''
	
	fields=sp_name.split()
	if len(fields) >= 2:
		sp_gen=fields[0]		
	else:
		sp_gen=sp_name
			
	return sp_gen


#*******************************************
# Fonction retournant le nom de l'espece (Genre espece)
#*******************************************
def getGenusSpecies(sp_name):

	'''Fonction retournant le nom de l'espece (Genre espece)	
	'''
	
	fields=sp_name.split()
	if len(fields) >= 2:
		sp_gen_sp=fields[0]+' '+fields[1]		
	else:
		sp_gen_sp=sp_name
			
	return sp_gen_sp


#*******************************************
# Fonction retournant le nom de l'espece (Genre espece)
#*******************************************
def getGenusSpeciesSP(sp_name):

	'''Fonction retournant le nom de l'espece (Genre espece et 3eme champ si SP.)	
	'''
	
	fields=sp_name.split()
	if len(fields) >= 2:
		if len(fields) >=3 and fields[1]=='SP.':
			sp_gen_sp=fields[0]+' '+fields[1]+' '+fields[2]			
		else:	
			sp_gen_sp=fields[0]+' '+fields[1]		
	else:
		sp_gen_sp=sp_name
			
	return sp_gen_sp


#*******************************************
# Fonction enlevant les doublons (ne garde qu'une seule souche) dans une liste de noms d'especes (
# 	Genre, Espece pour mode='genus_species' ou Genre si mode='genus'
# entree : liste de noms d'especes avec possibilite de redondance (souches)
# sortie : liste de noms d'especes non redondante => especes a utiliser pour la phylogenie
#	la 1ere des souches de l'espece rencontree dans le fichier est selectionnee arbitrairement
# sortie : liste des noms de souches non utilisees
#*******************************************
def deleteCopies_Sp(list_sp, mode='genus_species_sp'):
	'''Fonction enlevant les doublons (ne garde qu'une seule souche) dans une liste de noms d'especes 
	
	Enleve les doublons de meme Genre, Espece pour mode='genus_species' (i.e doublons de souches) ou Genre si mode='genus'
	entree : liste de noms d'especes avec possibilite de redondance (souches)
	sortie 1: liste de noms d'especes non redondante => especes a utiliser pour la phylogenie
	la 1ere des souches de l'espece rencontree dans le fichier est selectionnee arbitrairement
	sortie 2: liste des noms de souches non utilisees
	'''
	sp_nr=[]
	sp_notsel=[]	
    	sp_nr.insert(0, list_sp[0])
   
	i=1
	
	# Modes de selection : sur genus_species ou sur genus pour considerer seulement 1 seule espece par genre :
	if mode == 'genus_species':		
		while i < len(list_sp):	
        		s=list_sp[i]		
			code=0 # Aucun code n'a ete entre
        		pst=0							
        		for j in range(len(sp_nr)):
        	       	 	t=sp_nr[j]
				# Si l'espece est deja dans la liste new : test sur les 2 1ers champs (genre, espece) de s et t
				if (getGenusSpecies(s)==getGenusSpecies(t)): # Ne suffit pas !!! Il faut tester aussi sur code espece (cas de agrobacterium par ex)				
					pst+=1 # On prend note de la presence de l'espece dans la liste								
        		if (pst==0):
				sp_nr.append(s)
			else:
				sp_notsel.append(s)
		
			i+=1 # Fin while
		
	elif mode == 'genus':							
		while i < len(list_sp):	
        		s=list_sp[i]		
			code=0 # Aucun code n'a ete entre
        		pst=0							
        		for j in range(len(sp_nr)):
        	       	 	t=sp_nr[j]
				#gen_sp=getGenusSpecies(s)
				# Si l'espece est deja dans la liste new : test sur les 2 1ers champs (genre, espece) de s et t
				if (getGenus(s)==getGenus(t)): 				
					pst+=1 # On prend note de la presence de l'espece dans la liste		
        		if (pst==0):
				sp_nr.append(s)
			else:
				sp_notsel.append(s)
		
			i+=1 # Fin while
	
	elif mode == 'genus_species_sp':
		# Mode prenant en compte le fait qu'existe des especes "indeterminees" : 'Genre SP. bla' et 'Genre SP. blabla' seront considerees toutes les 2
		while i < len(list_sp):	
        		s=list_sp[i]		
			code=0 # Aucun code n'a ete entre
        		pst=0							
        		for j in range(len(sp_nr)):
        	       	 	t=sp_nr[j]
				# Si l'espece est deja dans la liste new : test sur les 2 1ers champs (genre, espece) de s et t
				if (getGenusSpeciesSP(s)==getGenusSpeciesSP(t)): 
					pst+=1 # On prend note de la presence de l'espece dans la liste								
        		if (pst==0):
				sp_nr.append(s)
			else:
				sp_notsel.append(s)
		
			i+=1 # Fin while
		
	
	
    	return (sp_nr, sp_notsel)#, dico_strains)

# Fin de la fonction deleteCopies_Sp()



#*******************************************
# Transforme une liste de noms de sequences (d'HOGENOM: "gene1_ESP1") en liste non redondante d'especes
# Renvoie la liste d'especes codees sans redondances
#*******************************************
def deleteCopies_SeqToSp(list_seq):#, dico_strains):
	'''Transforme une liste de noms de sequences (d'HOGENOM: "gene1_ESP1") en liste non redondante d'especes

	Renvoie la liste d'especes codees sans redondances
	'''


	list_sp=[]
	list_nr=[]
		
	for seq in list_seq:
		sp=seq.split(spsplitchar)[1]
		list_sp.append(sp)
	
	list_nr.append(list_sp[0])	
	i=1
	while i < len(list_sp):
        	s=list_sp[i]
        	pst=0
        	for j in range(len(list_nr)):
                	t=list_nr[j]
			# Si le code espece comporte 5 caracteres, il peut y avoir une difference sur le dernier pour la meme esp et souche (ex AGTU0 et AGTU1)
                	#if len(s)==5:
			#	if s[:-1]==t[:-1]:
			#		pst+=1
			#elif s==t: 
			#	pst+=1 # On prend note de la presence de l'element s dans la liste
			if s==t:
				pst+=1
		
        	if (pst==0):
			list_nr.append(s)
		i+=1
				
    	return (list_sp, list_nr)
	
# Fin de la fonction deleteCopies_SeqToSp(list_seq)

#~ #*******************************************
#~ # Fonction generant le fichier contenant les sequences en fasta pour le phylum et la famille donnee
#~ # (si le fichier fasta de la famille n'existe pas deja)
#~ # Retourne le nom du fichier fasta de la famille
#~ # REM : BD doit etre ouverte !!! 
#~ #*******************************************
#~ def getFastaSeqs(BD, wd, phylum, dico_phylums_species ,fam):
	#~ ''' Fonction generant le fichier contenant les sequences en fasta pour le phylum et la famille donnee

	#~ (si le fichier fasta de la famille n'existe pas deja)
	#~ Retourne le nom du fichier fasta de la famille
	#~ REM : BD doit etre ouverte !!! 
	#~ '''
	
	#~ #res=acnuc.acnucopen(BD)
	#~ #if res==0:
	
	#~ fasta_file=getPhylumFilename(wd, phylum, '_'+fam+'.fasta')
	#~ #if not os.path.exists(fasta_file):
	#~ if ((not os.path.exists(fasta_file)) or (os.path.getsize(fasta_file)==0) or (not os.path.exists(fasta_file+'-aln'))):	
		#~ sp_notsel=dico_phylums_species[phylum]['sp_notsel']
		#~ #print sp_notsel
		
		#~ dico_res={}
	
		#~ #print acnuc.alllistranks()
		#~ #print acnuc.getlistrank('fam_phyl')
		#~ #query='k='+fam+' et sp='+phylum+' et no sp=plasmid@'
		#~ query='k='+fam+' et sp='+phylum+' et o=nuclear' # Sequences plasmidiques mais pas organellaires prises en compte.
		
		
		#~ #print query
		#~ res_query=acnuc.proc_requete(query, 'fam_phyl')
		#~ #print acnuc.proc_requete(query, 'fam_phyl')
		#~ #print acnuc.proc_requete(query, fam)
				
		#~ # Les sequences des especes non choisies ne sont pas extraites :
		#~ for notsel in sp_notsel:
			#~ #print notsel
			#~ # Pour gestion des requetes trop longues en nb de caracteres : 
			#~ query='fam_phyl et no sp='+notsel
			#~ #query=fam+' et no sp='+notsel
			#~ #print query
			#~ res_query=acnuc.proc_requete(query, 'fam_phyl')
			#~ #res_query=acnuc.proc_requete(query, fam)
			#~ #print res_query
		
		#~ #print query
		#~ #res_query=acnuc.proc_requete(query, 'fam_phyl')
		#~ #res_query=acnuc.proc_requete(query, fam)
		#~ #print res_query
		
		#~ lrank=res_query['lrank'] 
			
		#~ #Generation du fichier fasta resultat de la requete :
		#~ open(fasta_file,'w').close() # Pour eviter d'etre en 'append'
		#~ acnuc.prep_extract('fasta', fasta_file, 'simple', lrank) # Attention en 'append' !!!!!
		#~ acnuc.extract_interrupt()
		#~ #acnuc.savelist(lrank, fam+'_test')
		#~ #print acnuc.alllistranks()
		
		#~ acnuc.releaselist(lrank)
		#~ #acnuc.acnucclose()
		
	#~ return fasta_file
	
	#~ #else:
	#~ #	print "probleme d'ouverture BD, code "+str(res)
	#~ #	return None
	
			

#~ # Fin de la fonction getFastaSeqs


#*******************************************
# Fonction prenant en argument un arbre au format phylip
# et retournant la liste des cles correspondantes.
#*******************************************
def extractSpeciesFromTree(tree_file):
	'''  Fonction prenant en argument un arbre au format phylip et retournant la liste des cles correspondantes.
	'''
	list_sp=[]

	#treefile=open('all_fam_ref_root.fasta.phylip_2_131fam_phyml_tree.txt', 'r')
	'''
	tree_file=open(tree_file, 'r')
	tree=tree_file.readlines()
	tree_file.close()
	'''
	tree=fileToLines(tree_file)
	
	tree=tree[0][:-2] # Recup du 1er arbre dans le fichier (sans '\n' et ';')


	tab=string.maketrans('','')
	tree=tree.translate(tab, '()') # Suppression des parentheses

	tree=tree.split(',')

	for x in tree:
		list_sp.append(x.split(':')[0])

	return list_sp	
	
#*******************************************
# Fonction prenant en argument un alignement au format phylip interleaved
# et retournant la liste des noms d'especes correspondantes. (code a 5 lettres)
#*******************************************	
def extractSpeciesFromPhylip(phy_file):	
	''' Fonction prenant en argument un alignement au format phylip interleaved et 
	retournant la liste des noms d'especes correspondantes. (code a 5 lettres)
	'''
	list_sp=[]	
	phy_file=open(phy_file, 'r')
	
	nb_sp=int(phy_file.readline().split()[0])

	i=nb_sp
	print nb_sp
	
	while i > 0:
		list_sp.append(phy_file.readline().split()[0])
		i-=1		
		
	return list_sp	


#*******************************************
# Fonction prenant en argument un fichier au format fasta
# et retournant la liste des noms d'especes correspondantes. (code a 5 lettres)
#*******************************************	
def extractSpeciesFromFasta(fasta_file):	
	'''Fonction prenant en argument un fichier au format fasta
	et retournant la liste des noms d'especes correspondantes. (code a 5 lettres)
	'''

	list_sp=[]
	
	fasta_file=open(fasta_file, 'r')
	
	parser=Fasta.RecordParser()
	iterator=Fasta.Iterator(fasta_file, parser)
	while 1:	
		record=iterator.next()
		if record is None:
			break
		sp=record.title.split()[0]
		if len(sp.split(spsplitchar)) > 1:
			sp=sp.split(spsplitchar)[1]
			list_sp.append(sp)
			
	fasta_file.close()
	#print list_sp
		
	return list_sp

#************************************************************
# Fonction realisant l'elagage d'un arbre de facon a obtenir 
# le meme echantillonage taxo que dans le fichier phy_file
# Prend en argument l'arbre et la liste des especes associeees
# et la liste des especes voulue en sortie (fam_sp_list)
#************************************************************
def pruneTree(wd, phylum, fam, phylum_sp_list, fam_sp_list, ref_tree):
	'''Fonction realisant l'elagage d'un arbre de facon a obtenir le meme echantillonage taxo que dans le fichier phy_file

	Prend en argument l'arbre et la liste des especes associeees et la liste des especes voulue en sortie (fam_sp_list)
	'''
	ref_tree_fam=getPhylumFilename(wd, phylum, '_'+fam+'_test_ref.tree')
	#print 'Arbre elague a creer : %s ' %(ref_tree_fam)
	
	# On supprime toutes les especes contenues dans la liste de la famille :
	for sp in fam_sp_list:
		# Cas ou fam_sp_list est inclu dans phylum_sp_list
		#if phylum_sp_list.__contains__(num):
		if phylum_sp_list.__contains__(sp):
			phylum_sp_list.remove(sp)
			
		# La famille presente des especes non contenues dans l'arbre de reference utilise : 
		else:
			# On doit supprimer les especes non presentes dans l'arbre de reference :
			# - dans le fichier phylip
			# - dans l'arbre de la famille
			print 'utilitaires.prunetree a finir...'
	
	# phylum_sp_list ne contient plus que les especes 'en trop', cad psentes dans l'arbre de ref et pas dans la famille
	shutil.copy(ref_tree, ref_tree_fam)
	
	for sp in phylum_sp_list:
		cmd='nopartial %s %s %s' %(ref_tree_fam, ref_tree_fam, sp)
		#print cmd
		os.system(cmd)
	
	return ref_tree_fam

# Fin de la fonction pruneTree


#*******************************************
# Enleve les doublons dans une liste 
# Renvoie la liste sans les redondances
#*******************************************
def deleteCopies(liste):
	'''Enleve les doublons dans une liste 	
	
	Renvoie la liste sans les redondances
	'''
	
	liste_nr=[]	
	liste_nr.append(liste[0])	
	i=1
	while i < len(liste):
        	s=liste[i]
        	pst=0
        	for j in range(len(liste_nr)):
                	t=liste_nr[j]
                	if (s==t): 
				pst+=1 # On prend note de la presence de l'element s dans la liste
		
        	if (pst==0):
			liste_nr.append(s)
		i+=1
				
    	return (liste_nr)
	


#********************************************************************************************
# Fonction recursive ecrivant toutes les combinaisons possibles de sequences a 1 seq/sp
# Au deb : dico avec str de donnees cf. loadDataHyp, la 1ere sp de la liste chainee, 
# une liste vide et un string vide
#********************************************************************************************
def comb_seq(dico, sp, hyp_list, str_comb):
	'''Fonction recursive ecrivant toutes les combinaisons possibles de sequences a 1 seq/sp

	Au deb : dico avec str de donnees cf. loadDataHyp, la 1ere sp de la liste chainee, une liste vide et un string vide
	'''

	if sp!=None:
		for seq in dico[sp]['seq']:
			str_comb_new=str_comb+' '+seq
			comb_seq(dico, dico[sp]['next'], hyp_list, str_comb_new)
	else:
		hyp_list.append(str_comb)
		#print hyp_list

#**************************************************************************************************************
# Fonction editant des fichiers au format fasta d'apres un dico comportant en cles le nom de sequences et en champ les sequences
# et d'apres une liste de noms de sequences choisies (sous-ensble de la liste des cles du dico)
#**************************************************************************************************************
def editFastaHypFiles(fasta_filename, dico, seq_list, hyp_num):
	'''Fonction editant des fichiers au format fasta d'apres un dico et une liste de noms de sequences
	
	Dico comportant en cles le nom de sequences et en champ les sequences
	et d'apres une liste de noms de sequences choisies (sous-ensble de la liste des cles du dico)
	'''

	hyp_filename=fasta_filename+'_hyp_'+str(hyp_num)+'.fasta'
	hyp_file=open(hyp_filename, 'w')	
	for seq_name in seq_list.split():
		hyp_file.write('> %s \n' %(seq_name) )
		hyp_file.write('%s\n' %(dico[seq_name]))
	
	hyp_file.close()
	
	return hyp_filename

#**************************************************************************************************************
# Fonction generant toutes les hypotheses d'orthologie possibles d'apres un fichier fasta apres alignement
#**************************************************************************************************************
def loadDataHyp(fasta_filename, fam_sp_list):

	'''Fonction generant toutes les hypotheses d'orthologie possibles d'apres un fichier fasta apres alignement
	'''

	fasta_file=open(fasta_filename, 'r')
	parser=Fasta.RecordParser()
	iterator=Fasta.Iterator(fasta_file, parser)	
	
	# Initialisation des str de donnees :
	dico_data={} # Structure de donnees (cf. liste chainee) pour stocker les sequences par especes
	dico_rec={} # Dico contenant les sequencs aux indices noms de sequences correspondant	
	i=0
	nb_sp=len(fam_sp_list)
	#print nb_sp
	for i in range(nb_sp):
		sp=fam_sp_list[i]
		dico_data[sp]={}
		dico_data[sp]['seq']=[]
		if i < (nb_sp-1):
			dico_data[sp]['next']=fam_sp_list[i+1]
		else:
			dico_data[sp]['next']=None
					
	# Remplissage des str de donnees :
	while 1:
		record=iterator.next()	
		if record is None:
			break
		
		seq_name=record.title.split()[0]
		sp=seq_name.split(spsplitchar)[1]
		dico_rec[seq_name]=mvWhiteSpaces(record.sequence)
		dico_data[sp]['seq'].append(seq_name)

	fasta_file.close()

	#print dico_rec
	#print dico_data.keys()
	hyp_list=[]
	str_comb=""
	#print fam_sp_list[0]
	comb_seq(dico_data, fam_sp_list[0], hyp_list, str_comb)
	
	nb_hyp=len(hyp_list) # Taille de hyp_list = nb d'hypotheses
	print "%d hypotheses" %nb_hyp

	i=0
	fasta_files_list=[]
	for i in range(nb_hyp):
		hyp=hyp_list[i]
		#print hyp
		fasta_files_list.append(editFastaHypFiles(fasta_filename, dico_rec, hyp, i+1))

	return fasta_files_list



#*******************************************
# Ne conserve que les elements presents dans les deux listes en argument 
# Renvoie la liste d'elements communs aux deux listes
#*******************************************
def commonLists(liste1, liste2):
	'''Conserve les elements communs a deux listes 	
	
	Renvoie la liste d'elements communs
	'''
	
	liste_com=[]
	if len(liste1) < len(liste2):	
		for el in liste1:
			if liste2.count(el)>0:
				liste_com.append(el)

	else:	
		for el in liste2:
			if liste1.count(el)>0:
				liste_com.append(el)
				
    	return liste_com
	

#******************************************* 
# Renvoie la liste d'elements presents dans une seule des listes
#*******************************************
def substractLists(liste1, liste2):
	'''Soustrait les elements de deux listes 
	
	Renvoie la liste d'elements uniques
	'''
	
	liste_sub=[]
	liste=liste1+liste2		
	for el in liste:
		if liste.count(el)==1:
			liste_sub.append(el)
					
    	return liste_sub




