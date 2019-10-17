#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""takes as input a fasta file with pantagruel genome code labels and turn it into an equivalent with organism names"""
import sys, os

fin = open(sys.argv[1], 'r')
fout = open(sys.argv[2], 'w')

dcodeorga = {}
with open(os.environ['database']+'/organism_codes.tab') as fcodeorga:
    for line in fcodeorga:
        lsp = line.rstrip('\n').split('\t')
        code = lsp[0]
        orga = lsp[1]
        stra = lsp[2]
        if not (stra in orga):
            orga += ' str. '+stra
        elif stra:
            if ' str. ' not in orga:
                orga = orga.replace(stra, 'str. '+stra)
        dcodeorga[code] = orga.replace('(', '').replace(')', '')

for line in fin:
    if line.startswith('>'):
        cdscode = line.strip('>\n')
        code = cdscode.split('_')[0]
        fout.write(">%s %s\n"%(cdscode, dcodeorga[code]))
    else:
        fout.write(line)

fin.close()
fout.close()
