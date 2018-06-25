#!/usr/bin/python

def coevol_score(mg_lineage_eventfreqs, dlineage_eventfreqs, nsample):
	cdef int f, f0, jointfreq
	jointfreq = 0
	for eid, f in mg_lineage_eventfreqs:
		f0 = dlineage_eventfreqs.get(eid)
		if f0: jointfreq += f0*f
	return float(jointfreq)/(nsample**2)

