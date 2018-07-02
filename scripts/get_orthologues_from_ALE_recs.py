#!/usr/bin/python

import sys, getopt, os
import tree2
from parseALErec import parseALERecFile

alemodel='dated'

for nfrec in lnfrec:
	spetree, subspetree, lrecgt, recgtlines, restrictlabs, dnodeevt = parseALERecFile(nfrec)
	
