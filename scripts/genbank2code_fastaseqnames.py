#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import sys, os
import multiprocessing, time

nflnfinfa = sys.argv[1]
nftranstab = sys.argv[2]
dirout = sys.argv[3]
if len(sys.argv)>4:
    nbcores = int(sys.argv[4])
else:
    nbcores = multiprocessing.cpu_count()

transnames = {}
with open(nftranstab, 'r') as ftranstab:
 for line in ftranstab:
  lsp = line.rstrip('\n').split()
  transnames[lsp[0]] = lsp[1]

with open(nflnfinfa, 'r') as flnfinfa:
  lnfinfa = [line.rstrip('\n') for line in flnfinfa]

def genbank2code(argtup):
	nfin, transnames, dirout, q = argtup
	nfout = os.path.join(dirout, os.path.basename(nfin).replace('.aln', '.codes.aln') )
	infile = open(nfin, 'r')
	outfile = open(nfout, 'w')
	for line in infile:
		if line.startswith('>'):
			lsp = line.strip('>\n').split()
			gbseqname = lsp[0].split('lcl|', 1)[-1]
			li = '>'+transnames[gbseqname]+'\n'
		else:
			li = line
		outfile.write(li)
	infile.close()
	outfile.close()
	q.put(nfout)
	return nfout

pool = multiprocessing.Pool(processes=nbcores)
manag = multiprocessing.Manager()
queue = manag.Queue(len(lnfinfa))

#~ result = pool.map_async(genbank2code, iter((nfinfa, transnames, dirout, queue) for nfinfa in lnfinfa))
pool.map(genbank2code, iter((nfinfa, transnames, dirout, queue) for nfinfa in lnfinfa))
#~ # monitor loop
#~ while not queue.full():
	#~ sys.stderr.write("\r%d\t/%d"%(queue.qsize(), len(lnfinfa)))
	#~ time.sleep(0.5)

#~ lnfoutfa = result.get()
#~ assert len(lnfoutfa) == len(lnfinfa)
