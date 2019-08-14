#!/usr/bin/python
import tree2

eventypes = 'DTSL'
outtaxlab = 'OUTGROUP'

def parseTERARec(frec, forcefreq=None, deadlab=outtaxlab):
	'''generator function that only parses the information related to the species tree
	
	NB: does not account for the change  of state (number of gene copies) over a species tree branch.
	For simplicity of representation, events are considered to happen at the begin (epper end) of the species tree branch.
	There is no implicit propagation of the gene copy number down the species tree; instead events of Speciations (S), 
	duplications (D) and incoming transfers (T) lead to an increment of the copy number in the recipient species,
	while losses (L) do not lead to decrement the count of gene copies (just no increment).
	'''
	dlevt = {et:[] for et in eventypes}
	dnodefreq = {}
	for line in frec:
		gfid, sevents, gchids = line.strip('\n').split(':')
		levents = sevents.split(';')
		for evtsup in levents:
			evt, sup = evtsup.split('@')
			freq = float(sup) if (not forcefreq) else float(forcefreq)
			evori, evtype, evdest1, evdest2 = (s.strip("'") for s in evt.split(','))
			if evtype=='S':
				dlevt['S'].append((evori, freq))
				dnodefreq[evdest1] = dnodefreq.setdefault(evdest1, 0.0) + 1.0
				dnodefreq[evdest2] = dnodefreq.setdefault(evdest2, 0.0) + 1.0
			if evtype=='D':
				dlevt['D'].append((evori, freq))
				dnodefreq[evdest1] = dnodefreq.setdefault(evdest1, 0.0) + 1.0
			elif evtype=='T':
				if evdest1==evori: evdon, evrec = (evdest1, evdest2)
				elif evdest2==evori: evdon, evrec = (evdest2, evdest1)
				dlevt['T'].append((evdon, evrec, freq))
			elif evtype=='SL':
				dlevt['S'].append((evori, freq))
				dlevt['L'].append((evdest1, freq))
				dnodefreq[evdest2] = dnodefreq.setdefault(evdest1, 0.0) + 1.0
			elif evtype=='TL':
				dlevt['T'].append((evori, evdest2, freq))
				dlevt['L'].append((evdest1, freq))
				dnodefreq[evdest2] = dnodefreq.setdefault(evdest2, 0.0) + 1.0
			elif evtype=='TTD':
				dlevt['T'].append((evori, deadlab, freq))
				dnodefreq[evdest2] = dnodefreq.setdefault(evdest1, 0.0) + 1.0
			elif evtype=='TFD':
				dlevt['T'].append((deadlab, evdest2, freq))
				dnodefreq[evdest1] = dnodefreq.setdefault(evdest2, 0.0) + 1.0
			elif evtype=='TLTD':
				dlevt['T'].append((evori, deadlab, freq))
				dlevt['L'].append((evdest1, freq))
			elif evtype=='TLFD':
				dlevt['T'].append((deadlab, evdest2, freq))
				# ignore losses in the dead/outgroup
				dnodefreq[evdest2] = dnodefreq.setdefault(evdest2, 0.0) + 1.0
	# only support single-reconciliation input files at the moment,
	# but return a sequence for consistency with analog functions
	return [(dnodefreq, dlevt)]	

def parseMowgliRec(frec):
	'''generator function that only parses the information related to the species tree
	
	take as input files containing (possibly single or multiple) recombinations, 
	separated by a line of dashes '--------' (8 or more of the character '-')
	'''
	def parseMowgliRecStreelines(frec):
		'''actual parser'''
		lsegments = []
		for line in frec:
			if line.startswith('--------') and lsegments:
			# delineate reconciliations by their '----' start line
				yield lsegments
				lsegments = []
			 
			if line.startswith('\t'):
				donrec, dates = line.strip('\t\n').split(' ')
				rec, don = donrec.strip('()').split(',')
				begend = [float(s) for s in dates.strip('[]').split(':')[0].split('-')]
				beg = begend[0]
				end = begend[-1] # this does not deal properly with dates in scientific notation
				lsegments.append((rec, don, beg, end))
		else:
			yield lsegments
	
	# interpret the information in terms of events
	for lsegments in parseTERARecFile(nfrec)
		dnodefreq = {}
		dlevt = {et:[] for et in eventypes}
		for seg in lsegments:
			rec, don, beg, end = seg
			# annotate gene occurence/frequency
			dnodefreq[rec] = dnodefreq.setdefault(rec, 0.0) + 1.0
			# determine the qgenerating event
			if not reftree[rec].is_child(reftree[don]):
				# transfer event
				dlevt['T'].append((don, rec, 1.0))
		yield (dnodefreq, dlevt)

def parseTERARecFile(nfrec, recformat='mowgli'):
	'''wrapper to parse reconciliations from ecceTERA, either in 'mowgli' (.mr) or 'tera' (.txt) format'''
	with open(nfrec, 'r') as frec:
		if recformat=='mowgli':
			return parseMowgliRec(frec)
		elif recformat=='tera':
			return parseTERARec(frec)
		else:
			raise ValueError, "wrong reconciliation format specification"
