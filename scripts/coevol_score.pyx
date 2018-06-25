#!/usr/bin/python

#~ def coevol_score(mg_lineage_eventfreqs, dlineage_eventfreqs, nsample):
	#~ cdef int f, f0, jointfreq
	#~ jointfreq = 0
	#~ for eid, f1 in mg_lineage_eventfreqs:
		#~ f0 = dlineage_eventfreqs.get(eid)
		#~ if f0: jointfreq += f0*f1
	#~ return float(jointfreq)/(nsample**2)


def coevol_score(currlineage_matches, nsample):
	cdef int rlocdsid, eid, f0, f1, jointfreq
	jointfreq = 0
	for rlocdsid, eid, f0, f1 in currlineage_matches:
		jointfreq += f0*f1
	return float(jointfreq)/(nsample**2)

