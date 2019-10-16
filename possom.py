# libraries
import os
import numpy as np
import pandas as pd
import ete3
import markov_clustering
import logging
import networkx
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# argument parser
arp = argparse.ArgumentParser()

# Add the arguments to the parser
arp.add_argument("-phy", "--phy",  required=True,  help="String. Folder with phylogenies")
arp.add_argument("-suf", "--suf",  required=True,  help="String. Suffix of phylogenies in folder")
arp.add_argument("-out", "--out",  required=True,  help="String. Prefix for output")
arp.add_argument("-ort", "--ort",  required=True,  help="String. Path to orthology file. Must be a two-column table (with tabs), with one gene per line: OG <tab> gene1")
arp.add_argument("-ani", "--ani",  required=True,  help="String. Which analysis to perform: \"main\" to analyse all genes, \"opti\" to find optimal inflation value")
arp.add_argument("-inf", "--inf",  required=False, default=1.1,      help="OPTIONAL: Floating. if analysis is \"main\", which inflation value to use? Default is 1.1")
arp.add_argument("-nopt", "--nopt", required=False, default=500,     help="OPTIONAL: Integer. if analysis is \"opti\", how many phylogenies should we examine for optimisation? Default is 500")
arp.add_argument("-print", "--print", required=False, default=False, help="OPTIONAL: Boolean (False/True). Print new tree with defined clusters?")
arp.add_argument("-split", "--split", required=False, default="_",   help="OPTIONAL: character to split species and sequence names. Default is \"_\", e.g. Human_genename. WARNING: use quotation marks, e.g. -split \"_\" or -split \"|\"")
arl = vars(arp.parse_args())

# input variables
phy_fo = arl["phy"]
phy_su = arl["suf"]
out_fn = arl["out"]
ort_fn = arl["ort"]
ani    = arl["ani"]
nopt   = int(arl["nopt"])
inf    = float(arl["inf"])
prb    = bool(arl["print"])
split_ch = arl["split"].replace("\"","")

# logging
logging.basicConfig(
	level=logging.DEBUG, 
	format="%(asctime)s [%(levelname)-5.5s]\t%(message)s",
	#handlers=[ logging.FileHandler("%s.log" % out_fn, mode="w"), logging.StreamHandler() ]
	handlers=[ logging.FileHandler("%s.log" % out_fn, mode="w") ]
	)

# prepare orthology from orthofinder
#os.system("awk '{ for (i=2; i <= NF; i++) { print $1\"\t\"$i  }}' %s | sed \"s/://\" > %s.orthology_preete.txt" % (ort_fn,out_fn))
ort = pd.read_csv("%s" % ort_fn, sep="\t", header=None)
if ort.shape[1] == 1: 
	ort["og"] = "NA"
ort.columns = ["node","og"]

# list of orthogroups
ort_lis = np.unique(ort["og"])



### FUNCTIONS ###

# loop through some orthogroups to find optimal inflation
def optimisation_loop(nopt=nopt):
	
	# loop through phylogenies
	n = 0
	inf_lis = np.zeros(nopt)
	mod_lis = np.zeros(nopt)

	ort_ran = np.random.choice(ort_lis, size=nopt)
	for n,phi in enumerate(ort_ran):

		# input name
		phy_fn = "%s/%s.%s" % (phy_fo,phi,phy_su)
		phy_id = phi.split(sep="/")[-1].split(sep=".")[0]

		# cluster phylogeny if you can, retrieve original clusters if you can't
		if os.path.exists(phy_fn):
			evd,_ = parse_phylo(phy_fn=phy_fn, phy_id=phy_id)
			inf_lis[n], mod_lis[n] = clusters_opt(phy_fn=phy_fn, phy_id=phy_id, evd=evd)
		else:
			inf_lis[n] = np.nan 
			mod_lis[n] = np.nan

		# counter
		if n % int(nopt/20) == 0 : print("#",n,"/",nopt)

	# calculate means
	print("#",nopt,"/",nopt)
	
	return inf_lis, mod_lis

# find optimal inflation for MCL
def optimise_inflation(matrix, start=1.1, end=2.5, step=0.1):
	I_lis = np.arange(start, end, step).tolist()
	Q_lis = np.zeros(shape=(len(I_lis),1))
	for n,I in enumerate(I_lis):
		result   = markov_clustering.run_mcl(matrix, inflation=I)
		clusters = markov_clustering.get_clusters(result)
		Q        = markov_clustering.modularity(matrix=result, clusters=clusters)
		Q_lis[n] = Q
	max_Q_index = np.argmax(Q_lis)
	# return inflation with maximum modularity and modularities array
	return I_lis[max_Q_index], Q_lis[max_Q_index]


# loop through all orthogroups and recluster them using 
# ETE + species overlap algorithm + MCL clustering
def orthogroup_loop():

	# loop through phylogenies
	n = 0
	l = len(ort_lis)
	for n,phi in enumerate(ort_lis):

		# input name
		phy_fn = "%s/%s.%s" % (phy_fo,phi,phy_su)
		phy_id = phi.split(sep="/")[-1].split(sep=".")[0]

		# cluster phylogeny if you can, retrieve original clusters if you can't
		if os.path.exists(phy_fn):
			evd,phy = parse_phylo(phy_fn=phy_fn, phy_id=phy_id)
			if len(evd) > 1:
				clu = clusters_mcl(phy_fn=phy_fn, phy_id=phy_id, evd=evd)
				if prb: print_tree(phy=phy, phy_fo=phy_fo, phi=phi, evc=clu)
			else:
				clu = clusters_nophylo(phy_id=phy_id, phy_fn=phy_fn)
		else:
			clu = clusters_nophylo(phy_id=phy_id, phy_fn=phy_fn)

		# create dataframe with old and new clusters, and all genes
		out = ort[ort["og"] == phy_id]
		out = pd.merge(out, clu, how="outer", on="node")
		out["og_cluster"] = out["og"] +"_"+ out["cluster"].astype(str)
		
		# save clusters
		if n == 0 : 
			out.to_csv("%s.orthology.csv" % out_fn, sep="\t", index=None, mode="w")
		if n > 0  : 
			out.to_csv("%s.orthology.csv" % out_fn, sep="\t", index=None, mode="a", header=False)
	
		# add counter
		if l > 20:
			if n % int(l/20) == 0 : print("#",n,"/",l)
		else:
			print("#",n,"/",l)
		

	# end triumphantly
	print("#",l,"/",l)


def print_tree(phy, phy_fo, phi, evc):

	for i in phy.get_leaves():
		c=evc[evc["node"] == i.name]["cluster"].values
		if c.size == 0: c="NA"
		else:           c=c[0]
		i.name = str(i.name) + "|cluster" + str(c) + "|"

	phy.write(outfile="%s/%s.ete_clusters.newick" % (phy_fo,phi))

# parse phylogenies with ETE to obtain a network-like table defining 
# orthologous relationships, using the species overlap algorithm
def parse_phylo(phy_fn, phy_id):

	# load input
	phy = ete3.PhyloTree("%s" % (phy_fn))
	logging.info("%s num nodes = %i" % (phy_id,len(phy)))
	# assign species names to tree
	phy.set_species_naming_function(lambda node: node.name.split(split_ch)[0] )
	# resolve polytomies in a random fashion
	phy.resolve_polytomy(recursive=True)
	# check if tree is rooted, apply midpoint root if unrooted
	phy_root = phy.get_tree_root()
	phy_outg = phy_root.get_children()
	is_root  = len(phy_outg) == 2
	if is_root:
		pass
		logging.info("%s Tree is rooted, pass" % phy_id)
	else: 
		logging.info("%s Tree is unrooted, apply midpoint root" % phy_id)
		phy_outgroup = phy.get_midpoint_outgroup()
		phy.set_outgroup(phy_outgroup)

	# find evolutionary events (duplications and speciations)
	evev = phy.get_descendant_evol_events(sos_thr=0)

	# create empty array for network edges
	evo    = np.empty((len(evev)*1000, 5), dtype="object")
	evo[:] = np.nan
	# loop through in and out seqs, create edge table with orthologous events
	n = 0
	for ev in evev:
		if ev.etype == "S":
			for ii in ev.in_seqs:
				for oi in ev.out_seqs:
					evo[n,0] = ii
					evo[n,1] = oi
					evo[n,2] = ev.branch_supports[0]
					evo[n,3] = ev.etype
					evo[n,4] = ev.sos
					n = n + 1

	evd = pd.DataFrame(evo).dropna()
	evd.columns = ["in_gene","out_gene","branch_support","ev_type","sos"]

	return evd, phy

# function to calculate optimal inflation for MCL for each phylogeny;
# i.e., that with the highest network modularity (Q)
def clusters_opt(phy_fn, phy_id, evd):
	
	# MCL clustering: create network
	logging.info("%s Create network" % phy_id)
	evd_e = evd[["in_gene","out_gene","branch_support"]] # define edges
	evd_n = networkx.convert_matrix.from_pandas_edgelist(evd_e, source="in_gene", target="out_gene", edge_attr="branch_support") # convert edges table into network
	evd_n_nodelist = [ node for i, node in enumerate(evd_n.node()) ] # list of nodes
	evd_m = networkx.to_scipy_sparse_matrix(evd_n, nodelist=evd_n_nodelist)
	# MCL clustering: find optimal inflation value
	inf,mod = optimise_inflation(matrix=evd_m) # find optimal inflation
	
	return inf, mod

# function to cluster a network-like table of orthologs (from ETE) using MCL
def clusters_mcl(phy_fn, phy_id, evd):

	# MCL clustering: create network
	logging.info("%s Create network" % phy_id)
	evou_e = evd[["in_gene","out_gene","branch_support"]]
	evou_n = networkx.convert_matrix.from_pandas_edgelist(evou_e, source="in_gene", target="out_gene", edge_attr="branch_support")
	evou_n_nodelist = [ node for i, node in enumerate(evou_n.node()) ]
	evou_m = networkx.to_scipy_sparse_matrix(evou_n, nodelist=evou_n_nodelist)
	# MCL clustering: run clustering
	# inf,_ = optimise_inflation(matrix=evou_m)
	logging.info("%s MCL clustering, inflation = %f" % (phy_id, inf))
	mcl_m  = markov_clustering.run_mcl(evou_m, inflation=inf)
	mcl_c  = markov_clustering.get_clusters(mcl_m)
	logging.info("%s MCL clustering, num clusters = %i" % (phy_id, len(mcl_c)))
	# markov_clustering.draw_graph(mcl_m, mcl_c, node_size=50, with_labels=True, edge_color="k", cmap="Accent")
	# MCL clustering: save output
	mcl_c_clu = [ i for i, cluster in enumerate(mcl_c) for node in cluster]
	mcl_c_noi = [ node for i, cluster in enumerate(mcl_c) for node in cluster]
	clu = pd.DataFrame( { 
		"node"    : [evou_n_nodelist[i] for i in mcl_c_noi],
		"cluster" : mcl_c_clu,
	}, columns=["node","cluster"])
	clu["cluster"] = clu["cluster"].astype(str)
	logging.info("%s MCL clustering, num clustered genes = %i" % (phy_id, len(clu)))

	return clu

# if the phylogeny can't be analysed with ETE (no phylogeny, not enough speciation events...), use
# the original orthogroups from orthofinder instead
def clusters_nophylo(phy_id, phy_fn):

	logging.info("%s Can't find phylogeny %s (small orthogroup?) OR can't run MCL (no speciations?), output original clusters instead" % (phy_fn,phy_id))
	clu = pd.DataFrame( { 
		"node"    : ort[ort["og"] == phy_id]["node"].values,
		"cluster" : np.nan
	}, columns=["node","cluster"])
	logging.info("%s No phylo/MCL, num clusters = %i" % (phy_id, 1))
	logging.info("%s No phylo/MCL, num clustered genes = %i" % (phy_id, len(clu)))

	return clu





### MAIN ####

# log input variables
logging.info("Input args: %r", arl)

# main analysis
if ani == "main":

	# run orthologs loop
	print("MODE: re-cluster orthogroups with ETE(SO)+MCL")
	logging.info("Analyse phylogenies in %s" % (phy_fo))
	orthogroup_loop()
	print("DONE")

elif ani == "opti":

	# run optimisation loop
	print("MODE: inflation optimisation with %i random phylogenies" % nopt)
	logging.info("Find optimal inflation")
	inf_lis, mod_lis = optimisation_loop(nopt=nopt)

	# plot optimisation histograms
	print("# median inflation is %f" % np.nanmedian(inf_lis))
	with PdfPages('%s.optimise_inflation.pdf' % out_fn) as pdf:
		# inflation
		plt.figure(figsize=(4,3))
		plt.title("Hist inflation I\nmedian = %f" % np.nanmedian(inf_lis))
		plt.hist(inf_lis[~np.isnan(inf_lis)])
		pdf.savefig(bbox_inches='tight')
		# modularity
		plt.figure(figsize=(4,3))
		plt.title("Hist modularity Q\nmedian = %f" % np.nanmedian(mod_lis))
		plt.hist(mod_lis[~np.isnan(mod_lis)])
		pdf.savefig(bbox_inches='tight')
		# close pdf
		plt.close()
	
	print("DONE")

else: 
	print("ERROR: Specify analysis with -ani (main or opti). See help with -h")

# end
logging.info("All done in %s" % (phy_fo))

