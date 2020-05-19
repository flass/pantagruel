#!/usr/bin/env Rscript --vanilla
library(ape)
cargs = commmadArgs()
for (carg in cargs){
	a = read.dna(carg, format='fasta')
	b = base.freq(a, freq=T)
#	names(b) = NULL
	print(b, quote=F)
}