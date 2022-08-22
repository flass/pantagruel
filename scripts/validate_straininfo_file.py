#!/usr/bin/env python2.7
import sys

forbiddenstuff = {'locus_tag_prefix':[None, '-_'], \
					'strain':[None, ' \/']}

nfin = sys.argv[1]
errprefix = "Error(s) in format of strain info file '%s': "%nfin
errsuffix = "Please edit the file accordingly."
dstrains = {}
expectedfields = set(['assembly_id','genus','species','strain','taxid','locus_tag_prefix'])
with open(nfin, 'r') as fin:
	headerline = fin.readline()
	header = headerline.rstrip('\n').split('\t')
	missingfields = expectedfields - set(header)
	extrafields = set(header) - expectedfields
	if missingfields:
		sextra = "\nExtra fields were detected: %s. Maybe the header was misspelled?"%repr(list(extrafields)) if extrafields else ''
		raise ValueError, "%s\nthe fields %s are missing from the header.%s\n%s"%(errprefix, repr(list(missingfields)), sextra, errsuffix)
	for fieldname in forbiddenstuff:
		fieldindex = header.index(fieldname)
		forbiddenstuff[fieldname][0] = fieldindex
	emptyfieldfound = 0
	charerrorfound = 0
	loctagprefixes = []
	for line in fin:
		lsp = line.rstrip('\n').split('\t')
		# check that every field has a value
		for fieldval, n in enumerate(lsp):
			if not fieldval:
				fieldname = header[n]
				emptyfieldfound += 1
				print "the field '%s' is empty for this entry:\n%s%s# # #\n"%(fieldname, headerline, line)
		for fieldname, fieldstuff in forbiddenstuff.iteritems():
			fieldindex, forbidenchars = fieldstuff
			fieldval = lsp[fieldindex]
			for fchar in forbidenchars:
				if (fchar in fieldval):
					print "the characters '%s' is forbidden in the '%s' field; rule broken at:\n%s%s# # #\n"%(forbidenchars, fieldname, headerline, line)
					charerrorfound += 1
			if fieldname == 'locus_tag_prefix':
				if fieldval in loctagprefixes:
					print "%s occurs more than once in the %s field"%(fieldval, fieldname)
				loctagprefixes.append(fieldval)
	
	errorlist = []
	if emptyfieldfound:
		errorlist.append("\t** empty fields were found")
	if charerrorfound:
		errorlist.append("\t** forbidden characters were found")
	if len(loctagprefixes) > len(set(loctagprefixes)):
		errorlist.append("\t** there are duplicate locus tag prefixe")
	if errorlist:
		raise ValueError, "%s\n%s\n%s"%(errprefix, '\n'.join(errorlist), errsuffix)
	