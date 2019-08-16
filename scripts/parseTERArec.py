#!/usr/bin/python
import tree2
import re

eventypes = 'DTSL'
outtaxlab = 'OUTGROUP'

realscinumpat = re.compile('[1-9]\.[0-9]{1,5}e-[0-9]{2}|0\.[0-9]{1,9}|0')

class TERARecSegment(object):
	"""store the string of reconciliation events carried by a reconciled gene tree branch, as modelled in TERA output file format"""
	
	def parseTERARecLine(recline, sgsep='_'):
		gfid, seventsgchids = recline.strip('\n').split(':', 1)
		if (sgsep in gfid): sevents = seventsgchids
		else: sevents, gchids = seventsgchids.split(':')
		levents = sevents.split(';')
		return (gfid, gchids, levents)
	
	def __init__(self, recline, sgsep='_'):
		parsedrecline = parseTERARecLine(recline, sgsep='_')
		self.gfid = parsedrecline[0]
		self.gchids = parsedrecline[1]
		self.sevents = parsedrecline[2]
		self.father = None # link to parent branch
		self.children = [] # links to child branches
		
	def set_father_segment(self, father):
		assert self.gfid in father.gchids
		self.father = father
		
	def set_child_segments(self, children):
		for c in children:
			assert c.gfid in self.gchids
			self.children.append(c)
		
	@staticmethod		
	def parseEventString(self):
		dlevt = {et:[] for et in eventypes}
		dnodefreq = {}
		levents = self.sevents.split(';')
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
			elif evtype.startswith('T'):
				if evdest1==evori: evdon, evrec = (evdest1, evdest2)
				elif evdest2==evori: evdon, evrec = (evdest2, evdest1)
				if evtype.endswith('TD'):
					evrec = deadlab
				elif evtype.endswith('FD'):
					evdon = deadlab
				# register event(s)
				if evtype.startswith('TL'):
					if evdon != deadlab:
						# loss event in the donor species; ignore in the dead/outgroup
						dlevt['L'].append((evdon, freq))
				elif evtype=='TLTD':
					dlevt['L'].append((evdon, freq))
				# transfer event
				dlevt['T'].append((evdon, evrec, freq))
				# register gene occurence
				if evrec != deadlab:
					# gain in the recipient species; ignore in the dead/outgroup
					dnodefreq[evrec] = dnodefreq.setdefault(evrec, 0.0) + 1.0
				if evdon != deadlab:
					# 'gain' in the donor species = maintenance of the lineage; ignore in the dead/outgroup
					dnodefreq[evdon] = dnodefreq.setdefault(evdon, 0.0) + 1.0
			elif evtype=='SL':
				dlevt['S'].append((evori, freq))
				dlevt['L'].append((evdest1, freq))
				dnodefreq[evdest2] = dnodefreq.setdefault(evdest1, 0.0) + 1.0
		return (dnodefreq, dlevt)

class TERAReconciliation(object):
	"""store the whole reconciliation of a gene tree, extracted from a TERA output file"""
		
	def __init__(self, recfile, sgsep='_'):
		self.dsegments = {}
		self.devents = {}
		self.dlevt = {}
		
		# parse the file, line per pile (each line a reconciliation segment corresponding to a reconciled gene tree branch)
		for line in recfile:
			rs = TERARecSegment(line, sgsep=sgsep)
			self.dsegments[rs.gfid] = rs
		# build the tree of segments = reconciliation
		for sid, segment in self.dsegments.items():
			segment.set_father_segment(self.dsegments[segment.gfid])
			segment.set_child_segments([self.dsegments[cid] for cid in segment.gchids])
	
			
	def shuntEventsToFromDead(self):

def parseTERARec(frec, forcefreq=None, deadlab=outtaxlab, sgsep='_'):
	"""generator function that only parses the information related to the species tree
	
	NB: does not account for the change  of state (number of gene copies) over a species tree branch.
	For simplicity of representation, events are considered to happen at the begin (epper end) of the species tree branch.
	There is no implicit propagation of the gene copy number down the species tree; instead events of Speciations (S), 
	duplications (D) and incoming transfers (T) lead to an increment of the copy number in the recipient species,
	while losses (L) do not lead to decrement the count of gene copies (just no increment).
	"""
	dlevt = {et:[] for et in eventypes}
	dnodefreq = {}
	for line in frec:
		gfid, seventsgchids = line.strip('\n').split(':', 1)
		if (sgsep in gfid): sevents = seventsgchids
		else: sevents, gchids = seventsgchids.split(':')
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
			elif evtype.startswith('T'):
				if evdest1==evori: evdon, evrec = (evdest1, evdest2)
				elif evdest2==evori: evdon, evrec = (evdest2, evdest1)
				if evtype.endswith('TD'):
					evrec = deadlab
				elif evtype.endswith('FD'):
					evdon = deadlab
				# register event(s)
				if evtype.startswith('TL'):
					if evdon != deadlab:
						# loss event in the donor species; ignore in the dead/outgroup
						dlevt['L'].append((evdon, freq))
				elif evtype=='TLTD':
					dlevt['L'].append((evdon, freq))
				# transfer event
				dlevt['T'].append((evdon, evrec, freq))
				# register gene occurence
				if evrec != deadlab:
					# gain in the recipient species; ignore in the dead/outgroup
					dnodefreq[evrec] = dnodefreq.setdefault(evrec, 0.0) + 1.0
				if evdon != deadlab:
					# 'gain' in the donor species = maintenance of the lineage; ignore in the dead/outgroup
					dnodefreq[evdon] = dnodefreq.setdefault(evdon, 0.0) + 1.0
			elif evtype=='SL':
				dlevt['S'].append((evori, freq))
				dlevt['L'].append((evdest1, freq))
				dnodefreq[evdest2] = dnodefreq.setdefault(evdest1, 0.0) + 1.0
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
				begend = [float(d) for d in reascinumpat.findall(dates)]
				beg = begend[0]
				end = begend[1] # ignore the possible 3rd member, always trailing after a colon and a repetion of the 1st and 2nd (when a duplication even?)
				if beg == end:
					# ignore this case as not sure what this [x,x:x] notation is about
					continue
				lsegments.append((rec, don, beg, end))
		else:
			yield lsegments
	
	# interpret the information in terms of events
	for lsegments in parseMowgliRecStreelines(nfrec):
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

def parseTERARecFile(nfrec, recformat='mowgli', **kw):
	'''wrapper to parse reconciliations from ecceTERA, either in 'mowgli' (.mr) or 'tera' (.txt) format'''
	with open(nfrec, 'r') as frec:
		if recformat=='mowgli':
			return parseMowgliRec(frec, **kw)
		elif recformat=='tera':
			return parseTERARec(frec, **kw)
		else:
			raise ValueError, "wrong reconciliation format specification"
