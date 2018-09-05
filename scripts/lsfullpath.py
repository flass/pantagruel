#!/usr/bin/python
import os, sys, glob

for target in sys.argv:
	if os.path.isdir(target):
		lfile = os.listdir(target)
		for f in lfile:
			print os.path.join(target, f)
	else:
		lpath = glob.glob(target)
			for f in lpath:
				print f
