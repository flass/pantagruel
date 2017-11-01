#!/usr/bin/python
import sys
from os import path
nfin = sys.argv[1]

dn = path.dirname(nfin)
bn = path.basename(nfin)
rad, ext = bn.rsplit('.', 1)
nfout = "%s/%s.nr.%s"%(dn, rad, ext)

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

print nfout
