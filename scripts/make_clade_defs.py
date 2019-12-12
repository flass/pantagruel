#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import tree2, sys
sys.setrecursionlimit(20000)

mincompare = 0
foutheader = ['', 'clade', 'sisterclade']

nfreftree = sys.argv[1]
if len(sys.argv)>2:
    mincompare = int(sys.argv[2])
    assert mincompare>0
    foutheader.append('backgroundclade')

reftree = tree2.read_check_newick(nfreftree)
nfout = nfreftree+"_clade_defs"
fout = open(nfout, 'w')
fout.write('\t'.join(foutheader)+'\n')
k = 0
for node in reftree:
    if len(node.children) != 2: continue
    for child in node.children:
        if child.nb_leaves() <= 1: continue
        focchildlab = child.label()
        if not focchildlab:
            focchildlab = "clade%d"%k
            k += 1
        child.edit_label(focchildlab)
        focchildleaflabset = set(child.get_leaf_labels())
        sischildleaflabset = set(child.go_brother().get_leaf_labels())
        outrow = [focchildlab, ','.join(sorted(focchildleaflabset)), ','.join(sorted(sischildleaflabset))]
        if mincompare>0:
            contrastnode = child.go_brother()
            contrastleaflabset = set(contrastnode.get_leaf_labels())
            while len(contrastleaflabset) < mincompare:
                if contrastnode.go_father():
                    contrastnode = contrastnode.go_father()
                    contrastleaflabset = set(contrastnode.get_leaf_labels()) - focchildleaflabset
                else:
                    break
            outrow.append(','.join(sorted(contrastleaflabset)))
                
        fout.write('\t'.join(outrow)+'\n')

fout.close()
reftree.write_newick(nfout+'.nwk', ignoreBS=True)