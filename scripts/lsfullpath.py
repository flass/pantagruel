#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import sys, glob

for target in sys.argv[1:]:
	lpath = glob.iglob(target)
	for f in lpath:
		print f
