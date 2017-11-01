#!/usr/bin/python
import os, sys
if len(sys.argv) > 1:
	target = sys.argv[1]
else:
	target = os.getcwd()
if len(sys.argv) > 2:
	pattern = sys.argv[2]
else:	
	pattern = ""
lfile = os.listdir(target)
for f in lfile:
	if not pattern:
		print "%s/%s"%(target, f)
	else:
		if pattern in f:
			print "%s/%s"%(target, f)
