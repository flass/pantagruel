#!/usr/bin/python

cdef double coevol_score(lineage_matches, unsigned long nsamplesq):
	cdef int f0, f1
	cdef unsigned long jf = 0
	cdef double jointfreq
	for tmatch in lineage_matches:
		f0 = tmatch[1]
		f1 = tmatch[2]
		jf += f0*f1
	jointfreq = jf
	return jointfreq/nsamplesq

#~ cdef chop_array(currlineage_matches, keepfields, chunksize=1000):
#~ 	cdef int tmatch_ff[2]
#~ 	cdef list lmatch_ffs = []
#~ 	for currlineage_match in currlineage_matches:
#~ 		# make an array of C-type int array
#~ 		tmatch_ff = currlineage_match[2:]
#~ 		lmatch_ffs.append(tmatch_ff)
#~ 	return lmatch_ffs

cpdef coevol_lineages(dbcur, lineage_id, nsamplesq, fetchsize=10000):
	"""generates a list of tuples containing pairs of lineage ids and the corresponding co-evolution score
	
	takes as input:
	- dbcur, a db cursor returning tuples from a querry with all integer values (rclocds_id, event_id, freq0, freq1), ordered by rclocds_id
	- lineage_id, the query lineage id
	- nsamplesq, the square of the number of sampled reconciliation, to scale the product of event frequencies into a probability
	- fetchsize, (optional) the number of rows fetched at once from the database
	"""
	if dbcur.rowcount == 0: return []
	cdef:
		unsigned long match_lineage_id, currlineage_id
		list  match_lineages, currlineage_matches = [], coevollineages = []
	match_lineages = dbcur.fetchmany(fetchsize)
	currlineage_id = match_lineages[0][0]
	while match_lineages:
		# seek boundaries of the lineage slices
		for k, tmatch_line_ff in enumerate(match_lineages):
			match_lineage_id = tmatch_line_ff[0]
			if match_lineage_id != currlineage_id:
				# finish parsing current lineage
				currlineage_matches += match_lineages[:k]
				# compute the lineages' co-evolution score
				coevollineages.append( (lineage_id, currlineage_id, coevol_score(currlineage_matches, nsamplesq)) )
				# next lineage
				currlineage_id = match_lineage_id
				currlineage_matches = []
				match_lineages = match_lineages[k:]
				break
		else:
			currlineage_matches += match_lineages
			match_lineages = dbcur.fetchmany(fetchsize)
	else:
		coevollineages.append( (lineage_id, currlineage_id, coevol_score(currlineage_matches, nsamplesq)) )
	return coevollineages
