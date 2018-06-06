#!/usr/bin/R
library(igraph)

#~ a = read.table('/home/flassall/PanteroDB_v0.3/07.compare_scenarios/bs_stem70_withinmedian35/replaceCCinGasinS-collapsePOPinSnotinG/ale_collapsed_undat/gene_event_matches/matching_events.tab')
a = read.table(commandArgs()[0])

q = quantile(a$V3, p=c(.99, .999, .9999, .99999, 1.0)
# top 1.0% have co-evol score in [11.0; 23.3] (2M pairwise associations)
# top 0.1% have co-evol score in [13.5; 23.3] (200k pairwise associations)
# top 0.01% have co-evol score in [15.7; 23.3] (20k pairwise associations)
# top 0.001% have co-evol score in [17.7; 23.3] (2k pairwise associations)
b = a[a$V3>=q['99.999%'],]
u1 = unique(b$V1)
u2 = unique(b$V2)
u = intersect(u1, u2)
# ~100 genes tested both ways, with strong association
g = graph_from_data_frame(b)

