#!/usr/bin/python2.7

import tree2, sys
sys.setrecursionlimit(10000)
nfcongt = sys.argv[1]
if nfcongt.endswith('.nwk'):
  congt = tree2.read_newick(nfcongt)
else:
  congt = tree2.read_nexus(nfcongt, returnDict=False, allLower=False)[0]
for node in congt:
  if not node.bs():
    node.set_bs(0.01)

congt.write_newick('%s.nwk'%nfcongt)
