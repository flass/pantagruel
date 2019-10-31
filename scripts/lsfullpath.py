#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import sys
import glob
import os

for target in sys.argv[1:]:
	if (('*' in target) or (('[' in target.replace('\[', '')) and (']' in target.replace('\]', '')))):
		lpath = glob.iglob(target)
	elif os.path.isdir(target):
		lpath = os.listdir(target)
	else:
		lpath = target
	for f in lpath:
		print f
