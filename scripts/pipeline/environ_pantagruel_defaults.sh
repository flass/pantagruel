#!/usr/bin/env bash
# Prokka annotation parameters (only relevant if custom genome assemblies are provided):
export assembler="somesoftware"
export seqcentre="somewhere"
export refgenus="Reference"
# species tree inference parameters
export ncorebootstrap=200
# gene tree inference parameters
export mainresulttag='rootedTree'
# gene trees collapsing DEFAULT values (used when -C option is NOT present in init call)
export cladesuppdef=70
export subcladesuppdef=35
export criteriondef='bs'
export withinfundef='median'
# bayesian gene tree estimation: separation of the many output files into separate folders
export ndiggrpfam=4
# gene tree/species tree reconciliation inference parameters
export ALEalgo='ALEml'
export ecceTERAalgo='amalgamate'
export recsamplesize=1000
# gene tree/species tree reconciliation parsing parameters for co-evolution analysis
export evtypeparse='ODTSL'
export minevfreqparse=0.1
export evtypematch='ODTS'
export minevfreqmatch=0.5
export minjointevfreqmatch=1.0
#~ export maxreftreeheight=0.25 # now a parameter with pantagruel -q

export ptgcitation="Lassalle F, Veber P, Jauneikaite E, Didelot X. Automated Reconstruction of All Gene Histories in Large Bacterial Pangenome Datasets and Search for Co-Evolved Gene Modules with Pantagruel.” bioRxiv 586495. doi: 10.1101/586495"
