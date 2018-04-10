#! /usr/local/bin/python
'''
Concatene des alignements quelque soit leur format et leur recouvrement taxonomique. 

La liste des fichiers d'alignements a concatener est passee en argument.
Author: Sophie Abby
'''

import lib_util
import utilitaires

import sys
import os
#import glob

GAP='-'
#GAP='_'
#aln_files=sys.argv[1]

def rebuildAln(sp_list_tot, aln_file):
	'''
	Fonction qui reconstitue un aln en integrant des especes pour lesquelles la proteine n'est pas forcement presente :
	ajout de gaps dans l'aln pour les especes en question, dont le nom est contenu dans sp_list. 
	'''
		
	aln=lib_util.AlnGenerator(aln_file)
	# 1_ Recuperation de la taille de l'aln 
	aln_lg=aln.get_aln_length()
	if aln_lg>0:
		# 2_ Pour chaque sp a ajouter : (add_sequence verifie d'abord si l'espece est dans l'alignement)
		for sp in sp_list_tot:		
			seq='%s'%GAP*aln_lg
			aln.add_sequence(sp,seq)
			#break
		aln.write_interleaved(aln_file+"_all_sp_new")
	
		return aln_file+"_all_sp_new"
	else:
		print "File %s NOT INCLUDED : no available aligned positions."%aln_file
		return None
		
def main():
	if len(sys.argv)<3:
		print "Usage : python concat.py liste_fichiers concatfilename."
		sys.exit( 1 )
	
	aln_files=sys.argv[1]	
	concatfile=sys.argv[2]
	if len(sys.argv)>3:
		nfrestrictsplist = sys.argv[3]
		with open(nfrestrictsplist, 'r') as frestrictsplist:
			lrestrictsp = [line.rstrip('\n') for line in frestrictsplist]
	else:
		lrestrictsp = []
		
	files=utilitaires.fileToLines(aln_files)
	#print files
	sp_list=[]
	aln_list=[]
	
	dico_aln={}
	for f in files:
		f=f[:-1]
		print f
		if os.path.exists(f):
			aln=lib_util.AlnGenerator(f)
			
			sp_list+=aln.get_species()
			aln_list.append(f)
			print "File %s included."%f	
		else:
			print "File %s NOT INCLUDED : does not exist."%f
	sp_list_nr=utilitaires.deleteCopies(sp_list)
	print "found %d species."%len(sp_list_nr)
	if lrestrictsp:
		print "restrict to species present in specified set of %d species"%len(lrestrictsp)
		sp_list_nr = list(set(sp_list_nr) & set(lrestrictsp))
		print "left %d species in the dataset"%len(sp_list_nr)
		
		##### does not work to restrict the alignement to species list specified in $3
	
	concat_list=[]	
	for aln in aln_list:
		res=rebuildAln(sp_list_nr, aln)
		#if res!=None:
		if res:
			concat_list.append(res)
	
	lib_util.concatAln(concat_list, concatfile, format="fasta")	
	
	return 0

main()
