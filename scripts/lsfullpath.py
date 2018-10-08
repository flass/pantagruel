#!/usr/bin/python
import sys, glob

for target in sys.argv[1:]:
	lpath = glob.iglob(target)
	for f in lpath:
		print f
