#!/usr/bin/python
import sys
from os import path
nfinaln = sys.argv[1]
nfoutaln = sys.argv[2]
nfrestrictspelist = sys.argv[3]
if len(sys.argv)>4:
	sep=sys.argv[4]
else:
	sep=None

with open(nfrestrictspelist, 'r') as frestrictspelist:
	lrestrictspe = [line.rstrip('\n') for line in frestrictspelist]

w = False
with open(nfinaln, 'r') as fin:
 with open(nfoutaln, 'w') as fout:
  for line in fin:
   if line.startswith('>'):
    lsp = line.strip('>\n').split()
    spe = lsp[0].split(sep)[0]
    if spe in lrestrictspe:
     w = True
    else:
     w = False
   if w: fout.write(line)
