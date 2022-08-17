#!/usr/bin/env python2.7
import sys

forbiddenstuff = {'locus_tag_prefix':[None, '-_'], \
                    'strain':[None, ' \/']}

nfin = sys.argv[1]
errprefix = "Error in format of strain info file '%s': "%nfin
errsuffix = "\nPlease edit the file accordingly."
dstrains = {}
expectedfields = set(['assembly_id','genus','species','strain','taxid','locus_tag_prefix'])
with open(nfin, 'r') as fin:
    header = fin.readline().rstrip('\n').split('\t')
    missingfields = expectedfields - set(header)
    extrafields = set(header) - expectedfields
    if missingfields:
        sextra = "\nExtra fields were detected: %s. Maybe the header was misspelled?"%repr(list(extrafields)) if extrafields else ''
        raise ValueError, "%s\nthe fields %s are missing from the header.%s%s"%(errprefix, repr(list(missingfields)), sextra, errsuffix)
    for fieldname in forbiddenstuff:
        fieldindex = header.index(fieldname)
        forbiddenstuff[fieldname][0] = fieldindex
    for line in fin:
        lsp = line.rstrip('\n').split('\t')
        for fieldname, fieldstuff in forbiddenstuff.iteritems():
            fieldindex, forbidenchars = fieldstuff
            fieldval = lsp[fieldindex]
            for fchar in forbidenchars:
                if (fchar in fieldval):
                    raise ValueError, "%s the characters '%s' is forbiden in the '%s' field. %s"%(errprefix, forbidenchars, fieldname, errsuffix)