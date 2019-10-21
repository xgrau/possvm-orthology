import argparse

# argument parser
arp = argparse.ArgumentParser()

# Add the arguments to the parser
arp.add_argument("-phy", "--phy",     required=True,                    help="String. Folder with phylogenies")
arp.add_argument("-out", "--out",     required=True,                    help="String. Prefix for output")
arp.add_argument("-mode", "--mode",   required=True,                    help="String. Which analysis to perform: \"single\" to analyse one gene tree (-phy); \"multi\" to analyse multiple trees in a folder (-phy and -suf required); \"opti\" to find optimal inflation value")
arp.add_argument("-pid", "--pid",     required=False, default="tree",   help="OPTIONAL: String. Name of gene family to analyse. REQUIRED if -mode single")
arp.add_argument("-suf", "--suf",     required=False, default="newick", help="OPTIONAL: String. Suffix of phylogenies in folder. REQUIRED if -mode multi")
arp.add_argument("-ort", "--ort",     required=False,                   help="OPTIONAL: String. Path to orthology file. REQUIRED if -mode multi. Must be a two-column table (with tabs), with one gene per line: OG <tab> gene1")
arp.add_argument("-root", "--root",   required=False, default="F",      help="OPTIONAL: Boolean (T/F). Are your trees rooted? If False, midpoint root is applied (default).")
arp.add_argument("-inf", "--inf",     required=False, default=1.1,      help="OPTIONAL: Floating. Which inflation value to use in MCL clustering? Default is 1.1")
arp.add_argument("-sos", "--sos",     required=False, default=0.0,      help="OPTIONAL: Floating. Which species overlap threshold used in ETE-SO? Default is 0.0")
arp.add_argument("-nopt", "--nopt",   required=False, default=500,      help="OPTIONAL: Integer. if analysis is \"opti\", how many phylogenies should we examine for optimisation? Default is 500")
arp.add_argument("-print", "--print", required=False, default="F",      help="OPTIONAL: Boolean (T/F). Print new tree with defined clusters?")
arp.add_argument("-split", "--split", required=False, default="_",      help="OPTIONAL: character to split species and sequence names. Default is \"_\", e.g. Human_genename. WARNING: use quotation marks, e.g. -split \"_\" or -split \"|\"")
arp.add_argument("-minbs", "--minbs", required=False, default=0,        help="OPTIONAL: Float. Minimum support for ortholog pairs. Orthologs linked by a tree branch with less support are dropped. Default is 0")
arp.add_argument("-refsp", "--refsp", required=False, default=None,     help="OPTIONAL: String. Indicate a reference species to use to annotate clusters of orthologs (adds an alphabetic label corresponding to a gene's cluster if that gene is orthologous to another from the ref sps)")
arl = vars(arp.parse_args())

# input variables
phy_fo   = arl["phy"]
out_fn   = arl["out"]
mod      = arl["mode"]
phy_id   = arl["pid"]
phy_su   = arl["suf"]
ort_fn   = arl["ort"]
nopt     = int(arl["nopt"])
inf      = float(arl["inf"])
sos      = float(arl["sos"])
do_print    = bool(arl["print"] == "T")
split_ch    = arl["split"].replace("\"","")
is_root     = bool(arl["root"] == "T")
min_support = arl["minbs"]
ref_sps     = arl["refsp"]

# libraries
import os
import numpy as np
import pandas as pd
import ete3
import markov_clustering
import logging
import networkx as nx
from networkx.algorithms import community
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import string


# os.chdir("/home/xavi/Documents/possom-orthology/")
# phy_fo = "/home/xavi/Documents/possom-orthology/test/single_genes/cyp_mosquitoes.newick"
# phy_fo = "/home/xavi/Documents/possom-orthology/test/anopheles_orthofinder/trees/OG0000000.newick"
# ref_sps="Anogam"
# phy_fo = "/home/xavi/Documents/possom-orthology/test/single_genes/adar_holozoa.newick"
# ref_sps="Hsap"
# out_fn = "ara"
# mod    = "single"
# phy_id = "adar"
# sos = 0
# inf = 1.1
# split_ch = "_"
# do_print = True
# nopt= 500
# phy_su="newick"
# ort_fn="res"
# is_root=False
# min_support=0


# logging
logging.basicConfig(
	level=logging.DEBUG, 
	format="%(asctime)s [%(levelname)-5.5s]\t%(message)s",
	handlers=[ logging.FileHandler("%s.log" % out_fn, mode="w"), logging.StreamHandler() ]
	#handlers=[ logging.FileHandler("%s.log" % out_fn, mode="w") ]
	)



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
			evs,_,_,_ = parse_phylo(phy_fn=phy_fn, phy_id=phy_id, is_root=is_root)
			inf_lis[n], mod_lis[n] = clusters_opt(phy_fn=phy_fn, phy_id=phy_id, evs=evs)
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


def print_tree(phy, out, evc, attributes, sep="|"):
	
	for i in phy.get_leaves():
		i_name = i.name
		for a in attributes:
			c=evc[evc["node"] == i_name][a].values
			if c.size == 0: c="NA"
			else:           c=c[0]
			i.name = str(i.name) + sep + a + str(c)
		i.name = str(i.name) + sep

	phy.write(outfile=out)

	# ts = ete3.TreeStyle()
	# ts.show_branch_support = True
	# ts.show_leaf_name = False
	# phy.render("%s.pdf" % out, tree_style=ts)


# parse phylogenies with ETE to obtain a network-like table defining 
# orthologous relationships, using the species overlap algorithm
def parse_phylo(phy_fn, phy_id, is_root):

	# load input
	phy = ete3.PhyloTree("%s" % (phy_fn))
	logging.info("%s num nodes = %i" % (phy_id,len(phy)))
	# assign species names to tree
	phy.set_species_naming_function(lambda node: node.name.split(split_ch)[0] )
	# resolve polytomies in a random fashion
	phy.resolve_polytomy(recursive=True)
	# check if tree is rooted, apply midpoint root if unrooted
	# NOT APPLIED: GLOBAL VARIBLE USED INSTEAD
	# phy_root = phy.get_tree_root()
	# phy_outg = phy_root.get_children()
	# is_root  = len(phy_outg) == 2
	if is_root:
		pass
		logging.info("%s Tree is rooted, pass" % phy_id)
	else: 
		logging.info("%s Tree is unrooted, apply midpoint root" % phy_id)
		phy_outgroup = phy.get_midpoint_outgroup()
		phy.set_outgroup(phy_outgroup)

	# ladderise phylogeny
	phy.ladderize()

	# list of genes in phylogeny
	phy_lis = phy.get_leaf_names()

	# find evolutionary events (duplications and speciations)
	evev = phy.get_descendant_evol_events(sos_thr=sos)

	# speciation events
	evs    = np.empty((len(evev)*len(evev), 5), dtype="object")
	evs[:] = np.nan
	n = 0
	for ev in evev:
		if ev.etype == "S":
			for ii in ev.in_seqs:
				for oi in ev.out_seqs:
					evs[n,0] = ii
					evs[n,1] = oi
					evs[n,2] = ev.branch_supports[0]
					evs[n,3] = ev.etype
					evs[n,4] = ev.sos
					n = n + 1
	evs = pd.DataFrame(evs).dropna()
	evs.columns = ["in_gene","out_gene","branch_support","ev_type","sos"]

	# duplications
	evd    = np.empty((len(evev)*len(evev), 5), dtype="object")
	# evd[:] = np.nan
	# n = 0
	# for ev in evev:
	# 	if ev.etype == "D":
	# 		for ii in ev.in_seqs:
	# 			for oi in ev.out_seqs:
	# 				evd[n,0] = ii
	# 				evd[n,1] = oi
	# 				evd[n,2] = ev.branch_supports[0]
	# 				evd[n,3] = ev.etype
	# 				evd[n,4] = ev.sos
	# 				n = n + 1
	# evd = pd.DataFrame(evd).dropna()
	# evd.columns = ["in_gene","out_gene","branch_support","ev_type","sos"]

	return evs, evd, phy, phy_lis


# function to calculate optimal inflation for MCL for each phylogeny;
# i.e., that with the highest network modularity (Q)
def clusters_opt(phy_fn, phy_id, evs):
	
	# MCL clustering: create network
	logging.info("%s Create network" % phy_id)
	evs_e = evs[["in_gene","out_gene","branch_support"]] # define edges
	evs_n = nx.convert_matrix.from_pandas_edgelist(evs_e, source="in_gene", target="out_gene", edge_attr="branch_support") # convert edges table into network
	evs_n_nodelist = [ node for i, node in enumerate(evs_n.node()) ] # list of nodes
	evs_m = nx.to_scipy_sparse_matrix(evs_n, nodelist=evs_n_nodelist)
	# MCL clustering: find optimal inflation value
	inf,mod = optimise_inflation(matrix=evs_m) # find optimal inflation
	
	return inf, mod


# function to cluster a network-like table of orthologs (from ETE) using MCL
def clusters_mcl(oid, evs, inf=inf):

	# MCL clustering: create network
	logging.info("%s Create network" % oid)
	evs_e = evs[["in_gene","out_gene","branch_support"]]
	evs_e = evs_e[evs_e["branch_support"] > min_support]
	evs_n = nx.convert_matrix.from_pandas_edgelist(evs_e, source="in_gene", target="out_gene", edge_attr="branch_support")
	evs_n_nodelist = [ node for i, node in enumerate(evs_n.node()) ]
	evs_m = nx.to_scipy_sparse_matrix(evs_n, nodelist=evs_n_nodelist)
	# MCL clustering: run clustering
	# inf,_ = optimise_inflation(matrix=evs_m)
	logging.info("%s MCL clustering, inflation = %f" % (oid, inf))
	mcl_m  = markov_clustering.run_mcl(evs_m, inflation=inf)
	mcl_c  = markov_clustering.get_clusters(mcl_m)
	logging.info("%s MCL clustering, num clusters = %i" % (oid, len(mcl_c)))
	# plt.figure(figsize=(10,10))
	#markov_clustering.draw_graph(mcl_m, mcl_c, node_size=20, with_labels=False, edge_color="silver", cmap="Accent")
	# plt.savefig('graph.pdf')
	# MCL clustering: save output
	mcl_c_clu = [ i for i, cluster in enumerate(mcl_c) for node in cluster]
	mcl_c_noi = [ node for i, cluster in enumerate(mcl_c) for node in cluster]
	clu = pd.DataFrame( { 
		"node"    : [evs_n_nodelist[i] for i in mcl_c_noi],
		"cluster" : mcl_c_clu,
	}, columns=["node","cluster"])
	clu["cluster"] = clu["cluster"].astype(str)
	logging.info("%s MCL clustering, num clustered genes = %i" % (oid, len(clu)))

	return clu


def ref_annot(clu, evs, ref_sps, syn_ref_nod):

	clu["sps"]     = clu["node"].apply(lambda c: c.split(split_ch)[0])
	evs["in_sps"]  = evs["in_gene"].apply(lambda c: c.split(split_ch)[0])
	evs["out_sps"] = evs["out_gene"].apply(lambda c: c.split(split_ch)[0])

	cluster_ref = []
	for n,noi in enumerate(clu["node"]):

		# which cluster?
		c = clu[clu["node"] == noi]["cluster"].values[0]

		# find reference sequences WITHIN THIS cluster 
		# and create a dictionary to alphabetic short codes
		ref_nodes = clu[(clu["sps"] == ref_sps) & (clu["cluster"] == c)]["node"].values
		ref_codes = list(string.ascii_letters[0:len(ref_nodes)])
		ref_dicti = dict()
		for m,r in enumerate(ref_nodes):
			ref_dicti[r] = ref_codes[m]

		# find if gene is orthologous to any ref sequences
		r1 = evs[(evs["in_gene"] == noi) & (evs["out_sps"] == ref_sps)]["out_gene"].values
		r2 = evs[(evs["out_gene"] == noi) & (evs["in_sps"] == ref_sps)]["in_gene"].values
		ra = np.unique(np.concatenate((r1,r2)))

		# if gene herself comes from a reference sequence, 
		# add it to the array of ref sequences
		if clu["sps"][n] == ref_sps:
			for syn in syn_ref_nod:
				if noi in syn:
					ra = np.unique(np.append(ra, list(syn)))

		# use cluster-specific alphabetic codes for ref sequences
		rc = np.sort([ ref_dicti[r] for r in ra ])
		rc = ''.join(rc)

		cluster_ref.append(c+rc)

	ix_ref_sps = np.where((clu["sps"] == ref_sps).values)
	cluster_ref_table = pd.DataFrame( { 
		"node" :        clu["node"].values[ix_ref_sps],
		"cluster" :     clu["cluster"].values[ix_ref_sps],
		"cluster_ref" : np.array(cluster_ref)[ix_ref_sps]
	} , columns=["node","cluster","cluster_ref"])
	
	return cluster_ref, cluster_ref_table


def find_monophyletic_refsps_expansion(phy, ref_sps):

	# list of species and gene names of descending leaves of each node in the tree
	node_species = phy.get_cached_content(store_attr="species")
	node_geneids = phy.get_cached_content(store_attr="name")

	# find sets of monophyletic sequences from the ref sps
	syn_seqs_refsps = list()
	for n,i in enumerate(node_species):
		if (len(node_species[i]) == 1) & (ref_sps in node_species[i]):
			syn_seqs_refsps.append(node_geneids[i])

	# remove redundancy in these sets (check for overlapping sets and merge into larger sets)
	syn_seqs_refsps_nr = syn_seqs_refsps
	for n,_ in enumerate(syn_seqs_refsps_nr):
		for m,_ in enumerate(syn_seqs_refsps_nr):
			if n != m:
				if len(syn_seqs_refsps_nr[n].intersection(syn_seqs_refsps_nr[m])) > 1:
					syn_seqs_refsps_nr[n] = syn_seqs_refsps_nr[n].union(syn_seqs_refsps_nr[m])
	syn_seqs_refsps_nr = np.unique(syn_seqs_refsps_nr)

	return syn_seqs_refsps_nr





# function to cluster a network-like table of orthologs (from ETE) using MCL
def clusters_mod(oid, evs, gene_list):

	# create network
	logging.info("%s Create network" % oid)
	evs_e = evs[["in_gene","out_gene","branch_support"]]
	evs_n = nx.convert_matrix.from_pandas_edgelist(evs_e, source="in_gene", target="out_gene", edge_attr="branch_support")
	# find communities using modularity
	logging.info("%s Find communities in network" % oid)
	evs_n_communities = list(community.greedy_modularity_communities(evs_n))
	#evs_n_communities = list(community.girvan_newman(evs_n))
	#evs_n_communities = list(community.asyn_lpa_communities(evs_n))
	#evs_n_communities = list(community.k_clique_communities(evs_n, 5))
	#evs_n_communities = list(community.asyn_fluidc(evs_n,2))
	clus_list    = np.zeros(len(gene_list))
	for n,noi in enumerate(gene_list):
		for com in range(len(evs_n_communities)):
			if noi in evs_n_communities[com]:
				clus_list[n] = int(com)+1
	# store clusters
	clu = pd.DataFrame( { 
		"node"    : gene_list,
		"cluster" : clus_list
	}, columns=["node","cluster"])
	clu["cluster"] = clu["cluster"].astype(int).astype(str)
	logging.info("%s Num clustered genes = %i" % (oid, len(clu)))
	logging.info("%s Num clusters = %i" % (oid, len(np.unique(clus_list))))

	return clu

# if the phylogeny can't be analysed with ETE (no phylogeny, not enough speciation events...), use
# the original orthogroups from orthofinder instead
def clusters_nomcl(oid, ort):

	logging.info("%s Can't run MCL (no speciations/no phylogeny), output original clusters instead" % (oid))
	clu = pd.DataFrame( { 
		"node"    : ort[ort["og"] == phy_id]["node"].values,
		"cluster" : np.nan
	}, columns=["node","cluster"])
	logging.info("%s No clustering, num clusters = %i" % (phy_id, 1))
	logging.info("%s No clustering, num genes = %i" % (phy_id, len(clu)))

	return clu





### MAIN ####

# log input variables
logging.info("Input args: %r", arl)

# main analysis

if mod == "single":

	print("MODE: find orthogroups in one phylogeny with ETE(SO)+MCL: %s" % phy_fo)
	logging.info("MODE: find orthogroups in one phylogeny with ETE(SO)+MCL: %s" % phy_fo)
	
	# read phylogeny, find speciation events, create network
	evs,evd,phy,phy_lis = parse_phylo(phy_fn=phy_fo, phy_id=phy_id, is_root=is_root)
	syn_ref_nod = find_monophyletic_refsps_expansion(phy=phy, ref_sps=ref_sps)

	# find clusters
	if len(evs) > 1:
		clu = clusters_mcl(oid=phy_id, evs=evs, inf=inf)
		#clu = clusters_mod(oid=phy_id, evs=evs, gene_list=phy_lis)
		if ref_sps is None:
			print_attributes = ["cluster"]
		else:
			clu["cluster_ref"], cluster_ref_table = ref_annot(clu=clu, evs=evs, ref_sps=ref_sps, syn_ref_nod=syn_ref_nod)
			cluster_ref_table.to_csv("%s.orthologs_refsps.csv" % out_fn, sep="\t", index=None, mode="w")	
			print_attributes = ["cluster","cluster_ref"]
		if do_print: 
			print_tree(phy=phy, out="%s.orthologs.newick" % out_fn, evc=clu, attributes=print_attributes, sep="|")
			os.system("nw_display %s.orthologs.newick -s -i visibility:hidden -d stroke:gray -b opacity:0 -w 1000 -v 12 > %s.orthologs.newick.svg" % (out_fn,out_fn))
			os.system("inkscape %s.orthologs.newick.svg --export-pdf=%s.orthologs.newick.pdf 2> /dev/null" % (out_fn,out_fn))
	else:
		print("No speciations in tree %s %s" % (phy_id, phy_fo))

	# save clusters
	clu.to_csv("%s.orthologs.csv" % out_fn, sep="\t", index=None, mode="w")
	
	# evolutionary events
	evs.to_csv("%s.orthologs_ete_speciation.csv" % out_fn, sep="\t", index=None, mode="w")
	# evd.to_csv("%s.orthologs_ete_duplication.csv" % out_fn, sep="\t", index=None, mode="w")



elif mod == "multi":

	# prepare orthology from orthofinder
	#os.system("awk '{ for (i=2; i <= NF; i++) { print $1\"\t\"$i  }}' %s | sed \"s/://\" > %s.orthology_preete.txt" % (ort_fn,out_fn))
	ort = pd.read_csv("%s" % ort_fn, sep="\t", header=None)
	if ort.shape[1] == 1: 
		ort["og"] = "NA"
	ort.columns = ["node","original_orthogroup"]

	# list of orthogroups
	ort_lis = np.unique(ort["original_orthogroup"])

	# run orthologs loop
	print("MODE: find orthogroups in tree collection with ETE(SO)+MCL: %s" % phy_fo)
	logging.info("MODE: find orthogroups in tree collection with ETE(SO)+MCL: %s" % phy_fo)

	# loop through phylogenies
	n = 0
	l = len(ort_lis)
	for n,phi in enumerate(ort_lis):

		# input name
		phy_fn = "%s/%s.%s" % (phy_fo,phi,phy_su)
		phy_id = phi.split(sep="/")[-1].split(sep=".")[0]

		# cluster phylogeny if you can, retrieve original clusters if you can't
		if os.path.exists(phy_fn):

			# read phylogeny, find speciation events, create network
			evs,evd,phy,phy_lis = parse_phylo(phy_fn=phy_fn, phy_id=phy_id, is_root=is_root)
			syn_ref_nod = find_monophyletic_refsps_expansion(phy=phy, ref_sps=ref_sps)

			# find clusters
			if len(evs) > 1:
				clu = clusters_mcl(oid=phy_id, evs=evs, inf=inf)
				if ref_sps is None:
					print_attributes = ["cluster"]
				else:
					clu["cluster_ref"], cluster_ref_table = ref_annot(clu=clu, evs=evs, ref_sps=ref_sps, syn_ref_nod=syn_ref_nod)
					cluster_ref_table.to_csv("%s.orthologs_refsps.csv" % out_fn, sep="\t", index=None, mode="w")	
					print_attributes = ["cluster","cluster_ref"]
				if do_print: 
					print_tree(phy=phy, out="%s.orthologs.newick" % out_fn, evc=clu, attributes=print_attributes, sep="|")
					os.system("nw_display %s.orthologs.newick -s -i visibility:hidden -d stroke:gray -b opacity:0 -w 1000 -v 12 > %s.orthologs.newick.svg" % (out_fn,out_fn))
					os.system("inkscape %s.orthologs.newick.svg --export-pdf=%s.orthologs.newick.pdf 2> /dev/null" % (out_fn,out_fn))
			else:
				clu = clusters_nomcl(oid=phy_id, ort=ort)
		else:
			clu = clusters_nomcl(oid=phy_id, ort=ort)

		# create dataframe with old and new clusters, and all genes
		out = ort[ort["original_orthogroup"] == phy_id]
		out = pd.merge(out, clu, how="outer", on="node")
		out["orthogroup"] = out["original_orthogroup"] +"_"+ out["cluster"].astype(str)
		
		# save clusters
		if n == 0 : 
			out.to_csv("%s.orthologs.csv" % out_fn, sep="\t", index=None, mode="w")
		if n > 0  : 
			out.to_csv("%s.orthologs.csv" % out_fn, sep="\t", index=None, mode="a", header=False)
	
		# add counter
		if l > 20:
			if n % int(l/20) == 0 : print("#",n,"/",l)
		else:
			print("#",n,"/",l)
		
	# end triumphantly
	print("#",l,"/",l)


elif mod == "opti":

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

