#!/usr/bin/python
import os, glob, sys, getopt
import tree2
import ptg_utils as ptg
from parseALErec import parseALERecFile, getOrthologues
import igraph
from itertools import combinations
import multiprocessing as mp
import cPickle as pickle

allmethods = ['strict', 'unicopy', 'mixed']

def summaryOGs(ogs, dlabs, N, verbose):
	n = len(ogs)
	if verbose: print "number of OGs:", n
	cov = sum([len(x) for x in ogs])
	if verbose: print "coverage of leaves:", cov, '/', N
	if cov != N:
		raise ValueError, "unclassified leaves:%s"%(repr(set(dlabs.values()) - set(reduce(lambda x,y: list(x)+list(y), ogs, []))))
	return str(n)

def writeRecGeneTreesColouredByOrthologs(lrecgt, logs, nfout, drevnexustrans, **kw):
	for recgt, ogs in zip(lrecgt, logs):
		ptg.colour_tree_with_leaf_groups(recgt, [[drevnexustrans[ll] for ll in og] for og in ogs])
		#~ ptg.colour_tree_with_constrained_clades(recgt, ogs)
	tree2.write_nexus(lrecgt, nfout=nfout, **kw)

def writeGraphCombinedOrthologs(nfoutrad, graph, clustering, recgt, llabs, drevnexustrans, colourCombinedTree, tag, **kw):
	# save the graph
	with open(nfoutrad+".orthologs.%s.pickle"%tag, 'w') as foutortgraphp:
		pickle.dump(graph, foutortgraphp, pickle.HIGHEST_PROTOCOL)
	# write out clustering = synthesis of OGs over the sample
	with open(nfoutrad+".orthologs.%s"%tag, 'w') as foutortgraphclust:
		foutortgraphclust.write('\n'.join(['%s\t%d'%(lab, mem) for lab, mem in zip(llabs, clustering.membership)])+'\n')
	# write a tree with coloured tips for graphic fun
	if colourCombinedTree:
		gcogs = [[graph.vs[v]['name'] for v in cluster] for cluster in clustering]
		writeRecGeneTreesColouredByOrthologs([recgt], [gcogs], nfoutrad+".orthologs.%s.nex"%tag, drevnexustrans, **kw)

def orthoFromSampleRecs(nfrec, outortdir, nsample=[], methods=['mixed'], \
                        foutdiffog=None, outputOGperSampledRecGT=True, colourTreePerSampledRecGT=False, \
                        graphCombine=None, majRuleCombine=None, **kw):
	""""""
	verbose = kw.get('verbose')
	fam = os.path.basename(nfrec).split('-', 1)[0]
	if verbose: print "\n# # # %s"%fam
	# collect the desired sample from the reconciliation file
	dparserec = parseALERecFile(nfrec, skipLines=True, skipEventFreq=True, nsample=nsample, returnDict=True)
	lrecgt = dparserec['lrecgt']
	if userefspetree:
		refspetree = dparserec['spetree']
	else:
		refspetree = None
	
	ddogs = {}
	dnexustrans = {}
	drevnexustrans = {}
	ltaxnexus = []
	llabs = []
	for i, recgenetree in enumerate(lrecgt):
		if nsample: g = nsample[i]
		else: g = i
		if verbose: print recgenetree
		if verbose: print "\n# # reconciliation sample %d"%g
		N = recgenetree.nb_leaves()
		dlabs = {}
		if set(['strict', 'mixed']) & set(methods):
			if verbose: print "\n# strict_ogs:\n"
			strict_ogs, unclassified, dlabs = getOrthologues(recgenetree, method='strict', refspetree=refspetree, dlabs=dlabs, **kw)
			n1 = summaryOGs(strict_ogs, dlabs, N, verbose)
		else:
			strict_ogs = unclassified = None; n1 = 'NA'
		if 'unicopy' in methods:
			if verbose: print "\n# unicopy_ogs:\n"
			unicopy_ogs, notrelevant, dlabs = getOrthologues(recgenetree, method='unicopy', refspetree=refspetree, dlabs=dlabs, **kw)
			n2 = summaryOGs(unicopy_ogs, dlabs, N, verbose)
		else:
			unicopy_ogs = None; n2 = 'NA'
		if 'mixed' in methods:
			if verbose: print "\n# mixed_ogs:\n"
			mixed_ogs, unclassified, dlabs = getOrthologues(recgenetree, method='mixed', strict_ogs=strict_ogs, unclassified=unclassified, refspetree=refspetree, dlabs=dlabs, **kw) #
			n3 = summaryOGs(mixed_ogs, dlabs, N, verbose)
		else:
			mixed_ogs = None; n3 = 'NA'
		
		if foutdiffog or verbose: 
			o12 = str(sum([int(o in strict_ogs) for o in unicopy_ogs])) if (strict_ogs and unicopy_ogs) else 'NA'
			o13 = str(sum([int(o in strict_ogs) for o in mixed_ogs])) if (strict_ogs and mixed_ogs) else 'NA'
			o23 = str(sum([int(o in unicopy_ogs) for o in mixed_ogs])) if (mixed_ogs and unicopy_ogs) else 'NA'
		if verbose:
			print "\n# summary:\n"
			print "overlap strict_ogs with unicopy_ogs:", o12
			print "overlap strict_ogs with mixed_ogs:", o13
			print "overlap unicopy_ogs with mixed_ogs:", o23
		if foutdiffog:
			foutdiffog.write('\t'.join([fam, str(g), n1, n2, n3, o12, o13, o23])+'\n')
		
		if colourTreePerSampledRecGT or colourCombinedTree:
			dnexuslabev2num = {}
			for l, leaf in enumerate(recgenetree.get_leaves()):
				sl = str(l)
				labev = leaf.label()
				lab = dlabs[labev]
				dnexuslabev2num[labev] = sl
				if i==0:
					dnexustrans[sl] = lab
					drevnexustrans[lab] = sl
					ltaxnexus = [dnexustrans[str(sl)] for sl in range(len(dnexustrans))]
			for node in recgenetree:
				if node.is_leaf():
					node.edit_label(dnexuslabev2num[node.label()])
				else:
					node.edit_label('')
		ddogs[g] = {'strict':strict_ogs, 'unicopy':unicopy_ogs, 'mixed':mixed_ogs}
		if verbose: print "\n# # # # # # # #"
		if i==0:
			# collect the leaf labels; just do once
			llabs = dlabs.values() 
			llabs.sort()
	
	R = len(lrecgt)
	gs = nsample if nsample else range(R)
	for method in methods:
		ltrees = []
		nfoutrad = os.path.join(outortdir, method, "%s_%s"%(fam, method))
		if colourTreePerSampledRecGT:
			logs = [ddogs[g][method] for g in gs]
			writeRecGeneTreesColouredByOrthologs(lrecgt, logs, nfoutrad+"_orthologous_groups.nex", drevnexustrans, \
				treenames=["tree_%d" for g in gs], ltax=ltaxnexus, dtranslate=dnexustrans, figtree=True)
		if outputOGperSampledRecGT:
			with open(nfoutrad+".orthologs.per_sampled_tree", 'w') as foutort:
				for g in gs:
					ogs = ddogs[g][method]
					foutort.write('\n'.join([' '.join(x) for x in ogs])+'\n#\n')
		
		if graphCombine or majRuleCombine:
			## for later output
			recgt0 = lrecgt[0] if colourCombinedTree else None 
			# could also use the ALE consensus tree, which has branch supports but has no lengths
			## first make a dict of edge frequencies
			dedgefreq = {}
			for g in gs:
				ogs = ddogs[g][method]
				for og in ogs:
					if len(og)==1:
						orfan = og[0] ; combo = (orfan, orfan)
						dedgefreq[combo] = dedgefreq.get(combo, 0) + 1
					else:
						# get all pairs of genes in the OG
						combogs = combinations(sorted(og), 2)
						# add the counts
						for combo in combogs:
							dedgefreq[combo] = dedgefreq.get(combo, 0) + 1
			## build a graph of connectivity of the genes in OGs, integrating over the sample
			gOG = igraph.Graph()
			gOG.add_vertices(len(llabs))
			gOG.vs['name'] = llabs
			if majRuleCombine:
				## make a majority rule unweighted graph
				mjgOG = gOG.copy()
				# drop edges with frequency below the threshold
				mjedges = []
				minfreq = majRuleCombine*R
				for e, f in dedgefreq.iteritems():
					# use strict majority (assuming the parameter majRuleCombine=0.5, the default) to avoid obtaining family-wide single components
					if f > minfreq: mjedges.append(e)
				# add the edges to the graph
				mjgOG.add_edges(mjedges)
				# find connected components (i.e. perform clustering)
				compsOGs = mjgOG.components()
				writeGraphCombinedOrthologs(nfoutrad, mjgOG, compsOGs, recgt0, llabs, drevnexustrans, colourCombinedTree, \
					tag="majrule_combined_%f"%majRuleCombine, ltax=ltaxnexus, dtranslate=dnexustrans, ltreenames=["tree_0"], figtree=True)
			if graphCombine:
				## make a full weighted graph
				# add the edges to the graph
				edges, freqs = zip(*dedgefreq.iteritems())
				gOG.add_edges(edges)
				gOG.es['weight'] = freqs
				# find communities (i.e. perform clustering)
				graphcommfun = getattr(igraph.Graph, "community_"+graphCombine)
				try:
					commsOGs = graphcommfun(gOG, weights='weight')
				except TypeError:
					try:
						commsOGs = graphcommfun(gOG, edge_weights='weight')
					except TypeError:
						commsOGs = graphcommfun(gOG)
				if isinstance(commsOGs, igraph.clustering.VertexDendrogram):
					commsOGs = commsOGs.as_clustering()
				assert isinstance(commsOGs, igraph.clustering.VertexClustering)
				writeGraphCombinedOrthologs(nfoutrad, gOG, commsOGs, recgt0, llabs, drevnexustrans, colourCombinedTree, \
					tag='graph_combined_%s'%graphCombine, ltax=ltaxnexus, dtranslate=dnexustrans, ltreenames=["tree_0"], figtree=True)

def usage(mini=False):
	basics = "Usage:\n python %s -i /path/to/input.reconciliation.dir -o /path/to/output.dir [OPTIONS]\n"%sys.argv[0]
	if mini: return basics
	s =  "Script for orthologous group (OG) assignment based on ALE reconciliation scenarios.\n\n"
	s += basics
	s += " Options:\n"
	s += "  --sample\t\t(comma-separated list of) index(es) of sampled reconciled gene trees to parse (default: all).\n"
	s += "  --methods\t\t(comma-separated list of) method(s) for orthologous group delineation."
	s += "\t\t\t\t\tCan be any non-empty combination of %s (default to 'mixed').\n"%repr(allmethods)
	s += "\t\t\t\t\t. NB: 'unicopy' method does not rely on ALE results other than the topology of sampled reconciled gene trees.\n"
	s += "  --use.species.tree\t\tadd species-tree based constraints on ortholo(only for methods 'unicopy' and 'mixed').\n"
	s += "  --threads\t\t\tnumber of parralel processes to run.\n"
	s += "  --verbose {0,1,2}\tverbose mode, from none to plenty.\n"
	s += "  -v\t\t\tequivalent to --verbose=1.\n"
	return s

if __name__=='__main__':
	
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:T:hv', ['input.reconciliation.dir=', 'output.dir=', 'ale.model=', 'sample=', \
														'methods=', 'max.frac.extra.spe=', 'use.species.tree', \
														'colour.sampled.trees', 'report.ogs.per.sampled.tree', 'summary.per.sampled.tree', \
														'majrule.combine=', 'graph.combine=', 'colour.combined.tree', \
														'threads=', 'verbose=', 'help'])
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	
	alerecdir = dopt.get('--input.reconciliation.dir', dopt.get('-i'))
	outortdir = dopt.get('--output.dir', dopt.get('-o'))
	if not (alerecdir and outortdir):
		print "Error: missing mandatory arguments"
		print usage(mini=True)
		sys.exit(2)
	
	userefspetree = dopt.get('--use.species.tree')
	verbose = int(dopt.get('--verbose', ('-v' in dopt)))
	ALEmodel = dopt.get('--ale.model', 'dated')
	summaryOverlap = ('--summary.per.sampled.tree' in dopt)
	outputOGperSampledRecGT = ('--report.ogs.per.sampled.tree' in dopt)
	colourTreePerSampledRecGT = ('--colour.sampled.trees' in dopt)
	colourCombinedTree = ('--colour.combined.tree' in dopt)
	trheshExtraSpe = float(dopt.get('--max.frac.extra.spe', 0.1))
	nsample = [int(k) for k in dopt.get('--sample', '').split(',') if k.isdigit()]
	methods = dopt.get('--methods', 'mixed').split(',')
	if not methods or set(methods) > set(allmethods):
		raise ValueError, "values for --methods must be any non-empty combination of %s (default to 'mixed').\n"%repr()
	graphCombine = dopt.get('graph.combine', 'fastgreedy')
	if not hasattr(igraph.Graph, "community_"+graphCombine):
		raise ValueError, "values for --graph.combine must be such that 'community_%val%' is a method of igraph.Graph objects, e.g. 'fastgreedy'; cf. http://igraph.org/python/doc/igraph.Graph-class.html"
	majRuleCombine = float(dopt.get('majrule.combine', 0.5))
	if not (majRuleCombine>0 and majRuleCombine<=1):
		raise ValueError, "values for --graph.combine must be a real within the interval ]0; 1]"
	nbthreads = int(dopt.get('--threads', dopt.get('-T', 1)))
	
	## main execution
	lnfrec = glob.glob('%s/*ale.ml_rec'%(alerecdir))

	for d in methods:
		pd = os.path.join(outortdir, d)
		if not os.path.isdir(pd):
			os.mkdir(pd)
	
	if summaryOverlap:
		if nbthreads>1:
			print "Warning: gathering per-tree summary of orthologous groups is DISABLED when running in parallel"
			foutdiffog = None	
		else:
			foutdiffog = open(os.path.join(outortdir, 'diff_ortho_methods_per_sampled_tree'), 'w')
			foutdiffog.write('\t'.join(['family', 'sampled_rec', 'nOG_strict', 'nOG_unicopy', 'nOG_mixed', 'overlap_strict_unicopy', 'overlap_strict_mixed', 'overlap_unicopy_mixed'])+'\n')
	else:
		foutdiffog = None
		
	if verbose and nbthreads>1:
		print "Warning: verbose mode is DISABLED when running in parallel"
		verbose = 0
			
	def orthoFromSampleRecsSetArgs(nfrec):
		orthoFromSampleRecs(nfrec, outortdir, nsample=nsample, ALEmodel=ALEmodel, \
							methods=methods, userefspetree=userefspetree, trheshExtraSpe=trheshExtraSpe, \
							graphCombine=graphCombine, majRuleCombine=majRuleCombine, colourCombinedTree=colourCombinedTree, \
							foutdiffog=foutdiffog, outputOGperSampledRecGT=outputOGperSampledRecGT, colourTreePerSampledRecGT=colourTreePerSampledRecGT, \
							verbose=verbose)
	
	# run sequentially
	if nbthreads==1:
		for nfrec in lnfrec:
			print nfrec
			orthoFromSampleRecsSetArgs(nfrec)
			if verbose: print " # # # # # # # \n"
		if foutdiffog: foutdiffog.close()
	# or run in parallel
	else:
		pool = mp.Pool(processes=nbthreads)
		res = pool.map(orthoFromSampleRecsSetArgs, lnfrec)
