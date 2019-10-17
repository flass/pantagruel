#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import sys
n=int(sys.argv[1])
N=int(sys.argv[2])
k = 0
while k*n < N:
  print '%d-%d'%((k*n)+1, min((k+1)*n, N)),
  k += 1
