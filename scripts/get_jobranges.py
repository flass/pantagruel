#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import sys
n = int(sys.argv[1])
N = int(sys.argv[2])
if len(sys.argv)>3:
	m = int(sys.argv[3])
else:
	m = 0

k = 0
while k*n < N:
	if (k==0 and m>0): a = m
	else: a = (k*n)+1
	b = min((k+1)*n, N)
	print '%d-%d'%(a, b),
	k += 1
