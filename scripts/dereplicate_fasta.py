#!/usr/bin/python
import sys
from os import path
nfin = sys.argv[1]
nfout = sys.argv[2]

sprot = set([])
w = False
with open(nfin, 'r') as fin:
 with open(nfout, 'w') as fout:
  for line in fin:
   if line.startswith('>'):
    lsp = line.strip('>\n').split()
    prot = lsp[0]
    if prot not in sprot:
     sprot.add(prot)
     w = True
    else:
     w = False
   if w: fout.write(line)
