#!/usr/bin/python
import sys, glob

for target in sys.argv:
	lpath = glob.iglob(target)
	for f in lpath:
		print f
