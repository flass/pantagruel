#!/bin/bash
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
# gene tree/species tree reconciliation inference parameters
export ALEalgo='ALEml'
export recsamplesize=1000
# gene tree/species tree reconciliation parsing parameters for co-evolution analysis
export evtypeparse='ST'
export minevfreqparse=0.1
export minevfreqmatch=0.5
export minjoinevfreqmatch=1.0
export maxreftreeheight=0.25
