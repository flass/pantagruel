#!/usr/bin/python

import gzip
import sys
try:
	fin = open(sys.argv[1], 'r')
except IOError, e:
	fin = gzip.open(sys.argv[1], 'rb')
	
fout = gzip.open(sys.argv[2], 'wb')

for line in fin:
	if line.startswith('>'):
		cdsname, product = line.strip('>\n').split(' ', 1)
		fout.write(">%s [locus_tag=%s] [protein=%s] [protein_id=%s] [gbkey=CDS]\n"%(cdsname, cdsname, product, cdsname))
	else:
		fout.write(line)
		
fin.close()
fout.close()
		
