#!/usr/bin/env Rscript
library(ape)
cargs = commandArgs(trailingOnly=T)
for (carg in cargs){
        a = read.dna(carg, format='fasta')
        b = base.freq(a, freq=T)
        cat(sprintf("%d %d %d %d\n", b['a'], b['c'], b['g'], b['t']))
}


