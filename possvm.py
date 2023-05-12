# libraries
import argparse
import sys
import os
os.environ['QT_QPA_PLATFORM']='offscreen'
import re
import numpy as np
import pandas as pd
import ete3
import logging

# argument parser
arp = argparse.ArgumentParser()

# Add the arguments to the parser
arp.add_argument("-i", "--in", required=True, help="Path to a phylogenetic tree in newick format. Each sequence in the tree must have a prefix indicating the species, separated from gene name with a split character. Default split character is \"_\", see --split for options.", type=str)
arp.add_argument("-o", "--out", required=False, default=None, help="OPTIONAL: String. Path to output folder. Defaults to same directory as input file.", type=str)
arp.add_argument("-p", "--phy",  required=False, default=None, help="OPTIONAL: String. Prefix for output files. Defaults to `basename` of input phylogeny. Default behaviour will never overwrite original files, because it adds suffixes.", type=str)
arp.add_argument("-r", "--ref",  required=False, default=None, help="OPTIONAL: String. Path to a table indicating reference gene names that can be used for orthogroup labeling. Format: geneid <tab> name.", type=str)
arp.add_argument("-refsps", "--refsps",  required=False, default=None, help="OPTIONAL: String. Comma-separated list of reference species that will be used for orthogroup labeling. If absent, all gene names present in the -r table will be considered.", type=str)
arp.add_argument("-s", "--sos", required=False, default=0.0, help="OPTIONAL: Float. Species overlap threshold used for orthology inference in ETE. Default is 0. Higher values (up to 1) result in more inclusive groupings.", type=float)
arp.add_argument("-method","--method", required=False, default="mcl", help="OPTIONAL: String. Clustering method. Options are `mcl` (MCL, default), `mclw` (MCL weighted by node supports), `louvain` (Louvain), or `lpa` (label propagation algorithm).", type=str)
arp.add_argument("-inflation","--inflation", required=False, default=1.5, help="OPTIONAL: Float. Inflation hyperparameter for MCL clustering. Only applicable if `method` is `mcl` or `mclw`. In practice, the most inflation-responsive method is `mclw`.", type=float)
arp.add_argument("-outgroup","--outgroup", required=False, default=None, help="OPTIONAL: String. Define a set of species that are treated as outgroups in the phylogeny, and excluded from orthology clustering. Can be a comma-separated list of species, or a file with one species per line. This option DOES NOT affect tree rooting, just orthology clustering. Disabled by default.", type=str)
arp.add_argument("-split", "--split", required=False, default="_", help="OPTIONAL: String to use as species prefix delimiter in gene ids, e.g. \"_\" for gene names formatted as speciesA_geneX. Defaults to \"_\".", type=str)
arp.add_argument("-itermidroot", "--itermidroot",  required=False, default=None, help="OPTIONAL: Integer. Turns on iterative midpoint rooting with INT iterations, which is used instead of the default midpoint rooting. Low numbers are recommended (e.g. 10 is often more than enough).", type=int)
arp.add_argument("-skiproot", "--skiproot",  required=False, action="store_false", help="OPTIONAL: Bool. Turns off tree rooting, in case your trees are already rooted.")
arp.add_argument("-skipprint","--skipprint", required=False, action="store_false", help="OPTIONAL: Bool. Turns off printing of annotated phylogeny in PDF format (annotated newick is still produced).")
arp.add_argument("-printallpairs","--printallpairs", required=False, action="store_true", help="OPTIONAL: Bool. Turns on the production of a table with pairwise orthology/paralogy relationships between all pairs of genes in the phylogeny (default behaviour is to only report pairs of orthologs).")
arp.add_argument("-min_support_node","--min_support_node", required=False, default=0, help="OPTIONAL: Float. Min node support to consider orthology relationships. If not set, all relationships are considered.", type=float)
arp.add_argument("-min_support_transfer","--min_support_transfer", required=False, default=None, help="OPTIONAL: Float. Min node support to allow transfer of labels from labelled to non-labelled groups in the same clade. If not set, this step is skipped.", type=float)
arp.add_argument("-clean_gene_names","--clean_gene_names", required=False, action="store_true", help="OPTIONAL: Bool. Will attempt to \"clean\" gene names from the reference table (see -r) used to create cluster names, to avoid very long strings in groups with many paralogs. Currently, it collapses number suffixes in gene names, and converts strings such as Hox2/Hox4 to Hox2-4. More complex substitutions are not supported.")
arp.add_argument("-cut_gene_names","--cut_gene_names", required=False, default=None, help="OPTIONAL: Integer. If set, will shorten cluster name strings to the given length in the PDF file, to avoid long strings in groups with many paralogs. Default is no shortening.", type=int)
arp.add_argument("-ogprefix","--ogprefix", required=False, default="OG", help="OPTIONAL: String. Prefix for ortholog clusters. Defaults to \"OG\".", type=str)
arp.add_argument("-spstree","--spstree", required=False, default=None, help="OPTIONAL: Path to a species tree. If this is provided, Possvm will use a species tree reconciliation algorithm instead of species overlap.", type=str)
arp.add_argument("-v","--version", action="version", version="%(prog)s 1.2")

arl = vars(arp.parse_args())


#########################
####### ARGUMENTS #######
#########################

# mandatory input variables
phy_fn = arl["in"]

# output folder
if arl["out"] is None:
	out_fn = os.path.dirname(phy_fn)
	if out_fn == "" :
		out_fn = "."
else:
	out_fn = arl["out"]

# check if out_fn exists, and create it if it doesn't
if not os.path.exists(out_fn):
    os.makedirs(out_fn)

# prefix for output files
if arl["phy"] is not None:
	phy_id = arl["phy"]
else:
	phy_id = os.path.basename(phy_fn)

# other parameters
sos                  = arl["sos"]
method               = arl["method"]
inflation            = arl["inflation"]
split_ch             = arl["split"].replace("\"","")
itermidroot          = arl["itermidroot"]
do_print             = arl["skipprint"]
do_allpairs          = arl["printallpairs"]
do_root              = arl["skiproot"]
min_support_transfer = arl["min_support_transfer"]
min_support_node     = arl["min_support_node"]
clean_gene_names     = arl["clean_gene_names"]
cut_gene_names       = arl["cut_gene_names"]
ogprefix             = arl["ogprefix"].replace("\"","")
spstree              = arl["spstree"]

# reference genes?
if arl["ref"] is not None:
	do_ref = True
	ref_fn = arl["ref"]
else:
	do_ref = False

# reference species?
if  arl["refsps"] is not None:
	refsps = arl["refsps"].split(",")
else:
	refsps = arl["refsps"]

# use an outgroup?
if arl["outgroup"] is not None:
	if os.path.exists(arl["outgroup"]):
		outgroup = pd.read_csv(arl["outgroup"], names=["species"])["species"].values
	else:
		outgroup = arl["outgroup"].replace(","," ").split()
else:
	outgroup = []


# select clustering method
valid_methods = ["mcl", "louvain", "lpa", "mclw"]
if method in valid_methods:
	clusters_function_string = "clusters_%s" % method
else:
	print("Error, invalid clustering method \'%s\'!" % method)
	print("Valid methods are: %s" % valid_methods)
	sys.exit()


# use species tree reconciliation?
if spstree is not None:
	do_sps_reconciliation = True
	phs = ete3.PhyloTree(spstree)
else:
	do_sps_reconciliation = False



#########################
####### FUNCTIONS #######
#########################

# logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)-5.5s]\t%(message)s", handlers=[ logging.StreamHandler() ] )

def write_tree(phy, out, evc, attributes, sep="|", do_print=True, cut_gene_names=120):
	
	logging.info("Print tree")
	
	phy_alter = phy.copy(method="newick-extended")
	
	for i in phy_alter.get_leaves():
		i_name = i.name
		for attribute in attributes:
			c=evc[evc["node"] == i_name][attribute].values
			if c.size == 0: 
				c="NA"
			else:  
				c=c[0]
			
			# cut if string is too long
			if cut_gene_names is not None:
				cut_gene_names = int(cut_gene_names)
				if len(c) > cut_gene_names+3:
					c=c[:cut_gene_names] + '...'
					
			i.name = str(i.name) + sep + str(c)
		i.name = str(i.name) + sep
	
	phy_alter.write(outfile=out)
	
	# print tree in pdf
	if do_print:
		ts = ete3.TreeStyle()
		ts.show_branch_support = True
		ts.show_leaf_name = False
		ts.complete_branch_lines_when_necessary = False
		ts.scale=120
		phy_alter.render("%s.pdf" % out, tree_style=ts)


# read in phylogeny and execute event parser to obtain table-like network of 
# orthologous relationships, using the species overlap algorithm
def parse_phylo(phy_fn, phy_id, do_root, do_allpairs, clusters_function_string, outgroup=outgroup, do_sps_reconciliation=do_sps_reconciliation):
	
	# load input
	phy = ete3.PhyloTree("%s" % (phy_fn))
	logging.info("%s num nodes = %i" % (phy_id,len(phy)))
	logging.info("%s clustering function is %s" % (phy_id,clusters_function_string))
	clusters_function = eval(clusters_function_string)
	if do_sps_reconciliation:
		logging.info("%s do species tree reconciliation" % (phy_id))
	
	# assign species names to tree
	phy.set_species_naming_function(lambda node: node.name.split(split_ch)[0] )
	
	# resolve polytomies (randomly)
	phy.resolve_polytomy(recursive=True)
	
	# try to find root if unrooted
	if do_root:
		
		# shall we do it with iterative midpoint rooting?
		if itermidroot is not None:
			
			niter = itermidroot
			num_evs_per_iter = np.zeros(niter)
			out_nod_per_iter = np.empty(niter,	dtype=object)
			
			# then, iterate to try to find better candidates
			phy_it = phy.copy(method="newick")
			phy_outgroup_it = phy_it.get_midpoint_outgroup()
			phy_it.set_outgroup(phy_outgroup_it)
			
			# select very short length to shorten every rooting candidate branch (based on in-tree branch distribution)
			dist_lengths = [ node.dist for node in phy.traverse("postorder") if node.dist > 0 ]
			shrunk_length = np.quantile(dist_lengths, 0.1)
			
			for roi in range(niter):
				
				# parse events and re-do clustering
				evs_it, _, _, phy_lis_it = parse_events(phy=phy_it, outgroup=outgroup, do_allpairs=False, min_support_node=min_support_node)
				clu_it = clusters_function(evs=evs_it, node_list=phy_lis_it, verbose=False)
				
				# store number of orthogroups in this particular iteration
				num_evs_per_iter[roi] = len(np.unique(clu_it["cluster"].values))
				out_nod_per_iter[roi] = phy_outgroup_it
				
				print("%s Iterative midpoint root | %i/%i | n OGs = %i" % (phy_id,roi+1,niter,num_evs_per_iter[roi]))
				# in subsequent iterations, shrink the previous outgroup branch, and try to find second-best candidate
				phy_outgroup_it.dist = shrunk_length
				phy_outgroup_it = phy_it.get_midpoint_outgroup()
				phy_it.set_outgroup(phy_outgroup_it)
			
			# select outgroup that minimises number of OGs (more agglomerative)
			phy_outgroup_ix = np.argmin(num_evs_per_iter)
			
			# outgroup node in iterated tree
			phy_outgroup_from_it = out_nod_per_iter[phy_outgroup_ix]
			phy_outgroup_descendants = [ t for t in phy_outgroup_from_it.get_leaf_names() ]
			
			# outgroup node in original tree
			phy_outgroup = phy.get_common_ancestor(phy_outgroup_descendants)
			
			if len(phy_outgroup_descendants) != len(phy_outgroup.get_leaf_names()):
				print("%s Iterative midpoint root found and impossible root, default to midpoint" % (phy_id))
				phy_outgroup_ix = 0
				phy_outgroup = phy.get_midpoint_outgroup()
				
			# set outgroup
			# print(phy_outgroup_ix, phy_outgroup)
			logging.info("%s Best root at iteration  | %i/%i | n OGs = %i" % (phy_id,phy_outgroup_ix+1,niter,num_evs_per_iter[phy_outgroup_ix]))
		
		# ...or shall we do it with simple midpoint rooting?
		else:
			
			# set outgroup using normal midpoint rooting
			logging.info("%s Midpoint root" % phy_id)
			phy_outgroup = phy.get_midpoint_outgroup()
			
		# set root
		phy.set_outgroup(phy_outgroup)
		
	# ignore rooting
	else: 
		
		pass
		logging.info("%s Skip rooting (assume tree is already rooted)" % phy_id)
		
	# ladderise phylogeny
	phy.ladderize()
	
	# parse events
	if do_sps_reconciliation:
		evs, eva, phy, phy_lis = parse_events_sps_reconciliation(phy=phy, phs=phs, outgroup=outgroup, do_allpairs=do_allpairs)
	else:
		evs, eva, phy, phy_lis = parse_events(phy=phy, outgroup=outgroup, do_allpairs=do_allpairs, min_support_node=min_support_node)
	clu = clusters_function(evs=evs, node_list=phy_lis)
	
	# output from event parsing
	return evs, eva, phy, phy_lis, clu


# parse phylogenies with ETE to obtain a network-like table defining 
# orthologous relationships, using the species overlap algorithm
def parse_events(phy, outgroup, do_allpairs, min_support_node=0):
	
	# list of genes in phylogeny
	phy_lis = phy.get_leaf_names()
	
	# find evolutionary events (duplications and speciations)
	evev = phy.get_descendant_evol_events(sos_thr=sos)
	
	# # speciation events
	# evs    = np.empty((len(evev)*len(evev), 5), dtype="object")
	# evs[:] = np.nan
	# n = 0
	# for ev in evev:
	# 	# retrieve list of sps in ingroup and outgroup:
	# 	sps_in = np.unique([ i.split(split_ch)[0] for i in ev.in_seqs ])
	# 	sps_ou = np.unique([ i.split(split_ch)[0] for i in ev.out_seqs ])
	# 	# check if node is a speciation node, or a duplication node where both descendant branches have exactly one species, and this is the same species
	# 	if (ev.etype == "S" or (ev.etype == "D" and len(sps_in) == len(sps_ou) == 1 and sps_in == sps_ou)) and ev.branch_supports[0] >= min_support_node:
	# 		for ii in ev.in_seqs:
	# 			for oi in ev.out_seqs:
	# 				evs[n,0] = ii
	# 				evs[n,1] = oi
	# 				evs[n,2] = ev.branch_supports[0]
	# 				evs[n,3] = ev.etype
	# 				evs[n,4] = ev.sos
	# 				n = n + 1
	# evs = pd.DataFrame(evs).dropna()
	# evs.columns = ["in_gene","out_gene","branch_support","ev_type","sos"]
	
	# speciation events
	# this new version avoids pre-creating a huge empty matrix that used to result in memory problems in large phylogenies
	def parse_ev(ev):
		sps_in = np.unique([ i.split(split_ch)[0] for i in ev.in_seqs ])
		sps_ou = np.unique([ i.split(split_ch)[0] for i in ev.out_seqs ])
		if (ev.etype == "S" or (ev.etype == "D" and len(sps_in) == len(sps_ou) == 1 and sps_in == sps_ou)) and ev.branch_supports[0] >= min_support_node:
			return(
				[{
					'in_gene':ii,
					'out_gene':oi,
					'branch_support':ev.branch_supports[0],
					'ev_type':ev.etype,
					'sos':ev.sos
				} for oi in ev.out_seqs for ii in ev.in_seqs]
			)
		else:
			return(None)
	outputs = [parse_ev(ev) for ev in evev]
	outputs = [x for x in outputs if x]
	evs = pd.DataFrame([item for sublist in outputs for item in sublist])
	
	# drop outgroup species, if any
	if len(outgroup) > 0:
		in_evs_sps =  [ i.split(split_ch)[0] in set(outgroup) for i in evs["in_gene"]  ]
		out_evs_sps = [ i.split(split_ch)[0] in set(outgroup) for i in evs["out_gene"] ]
		evs = evs.drop(np.where(np.logical_or(in_evs_sps, out_evs_sps))[0])
		phy_lis =  [ i for i in phy_lis if i.split(split_ch)[0] not in set(outgroup) ]
	
	# duplications and speciation events
	if do_allpairs:
		eva    = np.empty((len(evev)*len(evev), 5), dtype="object")
		eva[:] = np.nan
		n = 0
		for ev in evev:
			for ii in ev.in_seqs:
				for oi in ev.out_seqs:
					eva[n,0] = ii
					eva[n,1] = oi
					eva[n,2] = ev.branch_supports[0]
					eva[n,3] = ev.etype
					eva[n,4] = ev.sos
					n = n + 1
		eva = pd.DataFrame(eva).dropna()
		eva.columns = ["in_gene","out_gene","branch_support","ev_type","sos"]
	else:
		eva = []
	
	return evs, eva, phy, phy_lis


# parse phylogenies with ETE to obtain a network-like table defining 
# orthologous relationships, using the species overlap algorithm
def parse_events_sps_reconciliation(phy, phs, outgroup, do_allpairs):

	# list of genes in phylogeny
	phy_lis = phy.get_leaf_names()
	
	# find evolutionary events (duplications and speciations)
	recon_tree, evev = phy.reconcile(phs)
	
	# speciation events
	evs    = np.empty((len(evev)*len(evev), 5), dtype="object")
	evs[:] = np.nan
	n = 0
	for ev in evev:
		# retrieve list of sps in ingroup and outgroup:
		sps_in = np.unique([ i.split(split_ch)[0] for i in ev.in_seqs ])
		sps_ou = np.unique([ i.split(split_ch)[0] for i in ev.out_seqs ])
		# check if node is a speciation node, or a duplication node where both descendant branches have exactly one species, and this is the same species
		if (ev.etype == "S" or (ev.etype == "D" and len(sps_in) == len(sps_ou) == 1 and sps_in == sps_ou)):
			for ii in ev.in_seqs:
				for oi in ev.out_seqs:
					evs[n,0] = ii
					evs[n,1] = oi
					evs[n,2] = 0
					evs[n,3] = ev.etype
					evs[n,4] = 0
					n = n + 1
	evs = pd.DataFrame(evs).dropna()
	evs.columns = ["in_gene","out_gene","branch_support","ev_type","sos"]
	
	# drop outgroup species, if any
	if len(outgroup) > 0:
		in_evs_sps =  [ i.split(split_ch)[0] in set(outgroup) for i in evs["in_gene"]  ]
		out_evs_sps = [ i.split(split_ch)[0] in set(outgroup) for i in evs["out_gene"] ]
		evs = evs.drop(np.where(np.logical_or(in_evs_sps, out_evs_sps))[0])
		phy_lis =  [ i for i in phy_lis if i.split(split_ch)[0] not in set(outgroup) ]
	
	# duplications and speciation events
	if do_allpairs:
		eva    = np.empty((len(evev)*len(evev), 5), dtype="object")
		eva[:] = np.nan
		n = 0
		for ev in evev:
			for ii in ev.in_seqs:
				for oi in ev.out_seqs:
					eva[n,0] = ii
					eva[n,1] = oi
					eva[n,2] = ev.branch_supports[0]
					eva[n,3] = ev.etype
					eva[n,4] = ev.sos
					n = n + 1
		eva = pd.DataFrame(eva).dropna()
		eva.columns = ["in_gene","out_gene","branch_support","ev_type","sos"]
	else:
		eva = []
	
	return evs, eva, phy, phy_lis


# function to cluster a network-like table of orthologs (from ETE) using label propagation algorithm
def clusters_lpa(evs, node_list, verbose=True):

	import networkx as nx
	from networkx.algorithms import community
	
	if len(evs) > 0:
		# clustering: create network
		if verbose:
			logging.info("Create network")
		evs_e = evs[["in_gene","out_gene","branch_support"]]
		evs_n = nx.convert_matrix.from_pandas_edgelist(evs_e, source="in_gene", target="out_gene", edge_attr="branch_support")
		evs_n.add_nodes_from(node_list)
		
		# clustering: asynchronous label propagation
		if verbose:
			logging.info("Find communities LPA")
		clu_c = community.asyn_lpa_communities(evs_n, seed=11)
		clu_c = { frozenset(c) for c in clu_c }
		if verbose:
			logging.info("Find communities LPA num clusters = %i" % len(clu_c))
		clu_c_clu = [ i for i, cluster in enumerate(clu_c) for node in cluster ]
		clu_c_noi = [ node for i, cluster in enumerate(clu_c) for node in cluster ]
	
	else:
		
		if verbose:
			logging.info("There are no speciation events in this tree.")
		clu_c_noi = node_list
		clu_c_clu = [ i for i in range(len(node_list)) ]
	
	# clustering: save output
	clu = pd.DataFrame( { 
		"node"    : clu_c_noi,
		"cluster" : clu_c_clu,
	}, columns=["node","cluster"])
	if verbose:
		logging.info("Find communities LPA | num clustered genes = %i" % len(clu))
	
	return clu

# function to cluster a network-like table of orthologs (from ETE) using label propagation algorithm
def clusters_kclique(evs, node_list, verbose=True):
	
	import networkx as nx
	from networkx.algorithms import community
	
	if len(evs) > 0:
		# clustering: create network
		if verbose:
			logging.info("Create network")
		evs_e = evs[["in_gene","out_gene","branch_support"]]
		evs_n = nx.convert_matrix.from_pandas_edgelist(evs_e, source="in_gene", target="out_gene", edge_attr="branch_support")
		evs_n.add_nodes_from(node_list)
		
		# clustering: asynchronous label propagation
		if verbose:
			logging.info("Find communities k-clique")
		clu_c = community.k_clique_communities(evs_n, k=2)
		clu_c = { frozenset(c) for c in clu_c }
		if verbose:
			logging.info("Find communities k-clique num clusters = %i" % len(clu_c))
		clu_c_clu = [ i for i, cluster in enumerate(clu_c) for node in cluster ]
		clu_c_noi = [ node for i, cluster in enumerate(clu_c) for node in cluster ]
	
	else:
		
		if verbose:
			logging.info("There are no speciation events in this tree.")
		clu_c_noi = node_list
		clu_c_clu = [ i for i in range(len(node_list)) ]
	
	# clustering: save output
	clu = pd.DataFrame( { 
		"node"    : clu_c_noi,
		"cluster" : clu_c_clu,
	}, columns=["node","cluster"])
	if verbose:
		logging.info("Find communities k-clique | num clustered genes = %i" % len(clu))
	
	return clu

# function to cluster a network-like table of orthologs (from ETE) using Louvain
def clusters_louvain(evs, node_list, verbose=True):
	
	import networkx as nx
	
	if len(evs) > 0:
		
		# clustering: create network
		if verbose:
			logging.info("Create network")
		evs_e = evs[["in_gene","out_gene","branch_support"]]
		evs_n = nx.convert_matrix.from_pandas_edgelist(evs_e, source="in_gene", target="out_gene", edge_attr="branch_support")
		evs_n.add_nodes_from(node_list)
		
		# clustering: Louvain
		import community as community_louvain
		
		if verbose:
			logging.info("Find communities Louvain")
		clu_x = community_louvain.best_partition(evs_n)
		clu_c = {}
		for k, v in clu_x.items():
			clu_c[v] = clu_c.get(v, [])
			clu_c[v].append(k)
		if verbose:
			logging.info("Find communities Louvain | num clusters = %i" % len(clu_c))
		clu_c_noi = [ n for i,c in enumerate(clu_c) for n in clu_c[c] ]
		clu_c_clu = [ c for i,c in enumerate(clu_c) for n in clu_c[c] ]
		
	else:
		
		if verbose:
			logging.info("There are no speciation events in this tree.")
		clu_c_noi = node_list
		clu_c_clu = [ i for i in range(len(node_list)) ]
	
	# clustering: save output
	clu = pd.DataFrame( { 
		"node"    : clu_c_noi,
		"cluster" : clu_c_clu,
	}, columns=["node","cluster"])
	if verbose:
		logging.info("Find communities Louvain | num clustered genes = %i" % len(clu))
	
	return clu


# function to cluster a network-like table of orthologs (from ETE) using MCL
def clusters_mcl(evs, node_list, inf=inflation, verbose=True):
	
	import markov_clustering
	import networkx as nx
	
	if len(evs) > 0:
		
		# MCL clustering: create network
		if verbose:
			logging.info("Create network")
		evs_e = evs[["in_gene","out_gene","branch_support"]]
		evs_n = nx.convert_matrix.from_pandas_edgelist(evs_e, source="in_gene", target="out_gene", edge_attr="branch_support")
		evs_n.add_nodes_from(node_list)
		evs_n_nodelist = [ node for i, node in enumerate(evs_n.nodes()) ]
		evs_m = nx.to_scipy_sparse_array(evs_n, nodelist=evs_n_nodelist)
		# MCL clustering: run clustering
		if verbose:
			logging.info("MCL clustering, inflation = %.3f" % (inf))
		mcl_m  = markov_clustering.run_mcl(evs_m, inflation=inf)
		mcl_c  = markov_clustering.get_clusters(mcl_m)
		if verbose:
			logging.info("MCL clustering, num clusters = %i" % (len(mcl_c)))
		# MCL clustering: save output
		mcl_c_clu = [ i for i, cluster in enumerate(mcl_c) for node in cluster]
		mcl_c_noi = [ node for i, cluster in enumerate(mcl_c) for node in cluster]
		mcl_c_nod = [ evs_n_nodelist[i] for i in mcl_c_noi ]

	else:
		
		# MCL clustering: create network
		if verbose:
			logging.info("There are no speciation events in this tree.")
		mcl_c_nod = node_list
		mcl_c_clu = [ i for i in range(len(node_list)) ]
	
	# output
	clu = pd.DataFrame( { 
		"node"    : mcl_c_nod,
		"cluster" : mcl_c_clu,
	}, columns=["node","cluster"])
	if verbose:
		logging.info("MCL clustering, num clustered genes = %i" % (len(clu)))
	
	return clu



# function to cluster a network-like table of orthologs (from ETE) using MCL
def clusters_mclw(evs, node_list, inf=inflation, verbose=True):
	
	import markov_clustering
	import networkx as nx
	
	if len(evs) > 0:
		
		# MCL clustering: create network
		if verbose:
			logging.info("Create network")
		evs_e = evs[["in_gene","out_gene","branch_support"]]
		evs_e.columns = ["in_gene","out_gene","weight"]
		evs_n = nx.convert_matrix.from_pandas_edgelist(evs_e, source="in_gene", target="out_gene", edge_attr="weight")
		evs_n.add_nodes_from(node_list)
		evs_n_nodelist = [ node for i, node in enumerate(evs_n.nodes()) ]
		evs_m = nx.to_scipy_sparse_array(evs_n, nodelist=evs_n_nodelist)
		# MCL clustering: run clustering
		if verbose:
			logging.info("MCL weighted clustering, inflation = %.3f" % (inf))
		mcl_m  = markov_clustering.run_mcl(evs_m, inflation=inf) #why pruning threshold HAS to be zero?
		mcl_c  = markov_clustering.get_clusters(mcl_m)
		if verbose:
			logging.info("MCL weighted clustering, num clusters = %i" % (len(mcl_c)))
		# MCL clustering: save output
		mcl_c_clu = [ i for i, cluster in enumerate(mcl_c) for node in cluster]
		mcl_c_noi = [ node for i, cluster in enumerate(mcl_c) for node in cluster]
		mcl_c_nod = [ evs_n_nodelist[i] for i in mcl_c_noi ]
		
	else:
		
		# MCL clustering: create network
		if verbose:
			logging.info("There are no speciation events in this tree.")
		mcl_c_nod = node_list
		mcl_c_clu = [ i for i in range(len(node_list)) ]
	
	# output
	clu = pd.DataFrame( { 
		"node"    : mcl_c_nod,
		"cluster" : mcl_c_clu,
	}, columns=["node","cluster"])
	if verbose:
		logging.info("MCL clustering, num clustered genes = %i" % (len(clu)))
	
	return clu


# add a tag to cluster name (known genes within cluster)
def ref_tagcluster(clu, evs, ref, cluster_label="cluster", ref_spi=None, label_ref_node="node", label_if_no_annot="", clean_gene_names=False):
	
	logging.info("Add annotations: reference genes in each cluster")
	
	if ref_spi is None:
		ref_s = ref
	else:
		ref_sps = ref["gene"].apply(lambda c: c.split(split_ch)[0]).values
		ref_s = ref[np.isin(ref_sps, ref_spi)]
	
	cluster_tags = dict()
	cluster_list = np.unique(clu[cluster_label].values)
	for c in cluster_list:
		ref_is_node = np.isin(ref_s["gene"].values, clu[clu[cluster_label] == c][label_ref_node].values)
		ref_names = ref_s["name"].values [ ref_is_node ]
		cluster_tag = "/".join( natural_sort(np.unique([ str(i) for i in ref_names ])) )
		cluster_tags[c] = cluster_tag
	
	for c in cluster_list:
		if cluster_tags[c] == "":
			cluster_tags[c] = label_if_no_annot
	
	cluster_ref = [ cluster_tags[c] for c in clu[cluster_label].values ]
	
	if clean_gene_names:
		cluster_ref = [ sanitise_genename_string(r) for r in cluster_ref ]
	
	return cluster_ref	

# annotate orthology relationships to known genes
def ref_known_any(clu, evs, ref, syn_nod=None, label_if_no_annot="NA", label_ref_node="node", ref_spi=None):
	
	logging.info("Add annotations: orthology to known genes")
	
	if syn_nod is None:
		syn_nod = []
	
	if ref_spi is None:
		ref_s = ref
	else:
		ref_sps = ref["gene"].apply(lambda c: c.split(split_ch)[0]).values
		ref_s = ref[np.isin(ref_sps, ref_spi)]
	
	# find reference sequences 
	ref_genes = ref_s["gene"] [ np.isin(ref_s["gene"], clu["node"] ) ].values
	ref_names = ref_s["name"] [ np.isin(ref_s["gene"], clu["node"] ) ].values
	ref_nodes = clu["node"] [ np.isin( clu[label_ref_node], ref_genes ) ].values
	
	ref_dicti = dict()
	for m,r in enumerate(ref_genes):
		ref_dicti[r] = ref_names[m]
	
	node_ref = []
	node_sup = []
	for noi in clu["node"]:
		
		# find if gene is orthologous to any ref sequences
		r1 = evs[ ( evs["in_gene"] == noi  ) & np.isin( evs["out_gene"], ref_nodes ) ] ["out_gene"].values
		r2 = evs[ ( evs["out_gene"] == noi ) & np.isin( evs["in_gene"] , ref_nodes ) ] ["in_gene"].values
		ra = np.unique(np.concatenate((r1,r2)))
		ra_is_clustered = np.isin(element=ra, test_elements=ref_nodes)
		ra = ra[ra_is_clustered]
		
		# if gene is in reference, add it to ra as orthologous to itself
		if noi in ref_nodes:
			ra = np.append(ra, noi)
		
		# get reference name from dict
		rc,ixu = np.unique([ ref_dicti[r] for r in ra ], return_index=True)
		rc = '/'.join(rc)
		# get reference-node support from phylogeny
		rs = [ str(phy.get_common_ancestor([r, noi]).support) if r != noi else str(100.0) for r in ra ]
		rs = [ rs[ix] for ix in ixu ]
		rs = '/'.join(rs)
		
		# name of cluster
		node_ref.append(rc)
		node_sup.append(rs)
	
	# add NA string to empty fields
	node_ref = [ label_if_no_annot if c == "" else c for c in node_ref ]
	node_sup = [ label_if_no_annot if c == "" else c for c in node_sup ]
	
	return node_ref, node_sup


# natural alphanumeric sorting function
def natural_sort(l): 
	
	convert = lambda text: int(text) if text.isdigit() else text.lower() 
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
	return sorted(l, key = alphanum_key)


# retrieve statistical supports for each orthology cluster
def find_support_cluster(clu, phy, cluster_label="cluster"):
	
	# input and output lists
	cluster_list = np.unique(clu[cluster_label].values)
	cluster_support_dict = dict()
	
	# loop for each cluster
	for c in cluster_list:
		
		# nodes in cluster
		nodes_in_cluster = clu[clu[cluster_label] == c]["node"].values
		# find support in oldest node in cluster
		cluster_support_dict[c] = phy.get_common_ancestor(nodes_in_cluster.tolist()).support
	
	cluster_supports = [ cluster_support_dict[c] for c in clu[cluster_label].values ]
	
	return cluster_supports


# transfer cluster-level annotation to closely-related monophyletic groups
def find_close_monophyletic_clusters(clu, phy, ref_label="cluster_ref", ref_NA_label="NA", cluster_label="cluster", min_support_transfer=0, splitstring="/", exclude_level=None, exclude_label="NA"):
	
	logging.info("Add annotations: extend annotations to monophyletic groups")
	
	# input and output lists
	has_label_ix = np.where(clu[ref_label] != ref_NA_label)[0]
	has_label = clu["node"] [ has_label_ix ].values
	has_no_label_ix = np.where(clu[ref_label] == ref_NA_label)[0]
	clusters_labeled = np.unique(clu[cluster_label][has_label_ix].values)
	clusters_nonlabeled = np.unique(clu[cluster_label][has_no_label_ix].values)
	
	# exclude certain groups of sequences?
	if exclude_level is not None:
		clusters_nonlabeled = clu[cluster_label][has_no_label_ix].values[ clu[exclude_level][has_no_label_ix].values != exclude_label ]
	
	# extended annotations: init arrays with previous information
	extended_annots =   np.empty(shape=clu[ref_label].values.shape[0], dtype=object)
	extended_clusters = np.empty(shape=clu[ref_label].values.shape[0], dtype=object)
	
	num_extended = 0
	
	if len(clusters_labeled) > 0:
		
		# loop for each non-labeled cluster
		for ci in clusters_nonlabeled:
			
			# list of nodes in non-labeled cluster ci
			nodes_i = clu[clu[cluster_label] == ci]["node"].values.tolist()
			
			# if there is only one node, root is the same node (otherwise defaults to tree root!)
			if len(nodes_i) > 1:
				parent_i = phy.get_common_ancestor(nodes_i).get_ancestors()[0]
			elif len(nodes_i) == 1:
				parent_i = phy.get_leaves_by_name(nodes_i[0])[0].get_ancestors()[0]
			
			# check if descendants from parent group have references
			while True:
				
				# find leaves descending from the parent node
				parent_i_descendants = parent_i.get_leaf_names()
				parent_i_support = parent_i.support
				
				# if we have reached the root node, assume max support
				if parent_i.is_root():
					parent_i_support = 1e6
				
				# check if sister has refs
				has_refs = np.any( np.isin(element=parent_i_descendants, test_elements=has_label) )
				if has_refs and parent_i_support > min_support_transfer :
					
					descendants_with_ref_ix = np.where( np.isin(element=clu["node"].values, test_elements=parent_i_descendants) )[0]
					nodes_in_descendants    = clu["node"] [ np.intersect1d ( descendants_with_ref_ix, has_label_ix ) ].values
					clusters_in_descendants = clu[cluster_label] [ np.intersect1d ( descendants_with_ref_ix, has_label_ix ) ].values
					annots_in_descendants   = clu[ref_label] [ np.intersect1d ( descendants_with_ref_ix, has_label_ix ) ].values
					nodes_in_descendants    = [ node for node in nodes_in_descendants ] # "flatten" list because it's broken for some reason
					clusters_in_descendants_list = np.unique(clusters_in_descendants)
					annots_in_descendants_list = np.unique(annots_in_descendants)
					num_extended = num_extended + 1
					
					break
				
				# if not, visit upper node, and go back to checking out its descendants
				else:
					parent_i = parent_i.get_ancestors()[0]
				
			# where to assign new labels?
			needs_new_label_ix = np.where( np.isin(element=clu["node"].values, test_elements=nodes_i) )[0]
			
			# reorder labels (alphabetically)
			flatten = lambda l: [item for sublist in l for item in sublist]
			annots_in_descendants_list = [t.split(splitstring) for t in annots_in_descendants_list]
			annots_in_descendants_list = np.unique(natural_sort(flatten(annots_in_descendants_list)))
			
			# join into a single string:
			clusters_in_descendants_string = "/".join( natural_sort([ str(i) for i in clusters_in_descendants_list ]) )
			annots_in_descendants_string = "/".join( natural_sort(  [ str(i) for i in annots_in_descendants_list   ]) )
			
			# store
			extended_clusters [ needs_new_label_ix ] = clusters_in_descendants_string
			extended_annots [ needs_new_label_ix ]   = annots_in_descendants_string
	
	logging.info("Add annotations: extend annotations to monophyletic groups | %i labels transferred" % num_extended)
	
	return extended_clusters, extended_annots


# attempt to shorten annotations (will collapse similarly named annotations)
def sanitise_genename_string(string, splitstring="/"):

	# gent gene names
	names = sorted(string.split(splitstring))
	# find prefixes (non-numeric characters at the beginning of gene name)
	prefixes = [ re.findall(r'^[^\d]+', name) or [""] for name in names ]
	# find gene number suffixes (numeric characters at the end of gene name)
	numbers =  [ re.findall(r'\d+[A-Za-z]*$', name) or [""] for name in names ]

	# flatten lists
	flatten = lambda l: [item for sublist in l for item in sublist]
	numbers = flatten(numbers)
	prefixes = flatten(prefixes)

	new_names = np.empty(shape=len(names), dtype=object)

	for n,prefix in enumerate(prefixes):
		# get numbers from identical prefixes
		same_prefixes_ix = np.where(np.isin(element=prefixes, test_elements=prefix))[0]
		same_prefixes_numbers = np.array(numbers)[same_prefixes_ix]
		same_prefixes_numbers_clean = natural_sort([ n for n in same_prefixes_numbers if n ])
		suffix_string = "-".join(same_prefixes_numbers_clean)
		new_names[n] = "".join([prefix, suffix_string])

	new_names = np.unique(new_names)
	new_names = "/".join(new_names)

	return new_names


# classify pairwise gene relationships according to paralogy/orthology, same/different species, and same/different cluster
def annotate_event_type(eva, clu, clutag="cluster_name", split_ch=split_ch):
	
	eva = pd.merge(left=eva, right=clu, left_on="in_gene", right_on="node", how="left")
	eva = pd.merge(left=eva, right=clu, left_on="out_gene", right_on="node", how="left")
	eva["in_sps"] = eva["in_gene"].apply(lambda c: c.split(split_ch)[0]).values
	eva["out_sps"] = eva["out_gene"].apply(lambda c: c.split(split_ch)[0]).values
	
	# empty
	evtype = np.empty(len(eva), dtype=object)
	# orthologs
	evtype[ ( eva["ev_type"] == "S" ) & ( eva["cluster_name_x"] == eva["cluster_name_y"] )& ( eva["in_sps"] != eva["out_sps"] ) ] = "ortholog_int"
	evtype[ ( eva["ev_type"] == "S" ) & ( eva["cluster_name_x"] != eva["cluster_name_y"] )& ( eva["in_sps"] != eva["out_sps"] ) ] = "ortholog_ext"
	# outparalog_ext (from different species and an external orthology group)
	evtype[ ( eva["ev_type"] == "D" ) & ( eva["cluster_name_x"] == eva["cluster_name_y"] ) & ( eva["in_sps"] != eva["out_sps"] ) ] = "outparalog_int"
	# outparalog_int (from different species and the same orthology group)
	evtype[ ( eva["ev_type"] == "D" ) & ( eva["cluster_name_x"] != eva["cluster_name_y"] ) & ( eva["in_sps"] != eva["out_sps"] ) ] = "outparalog_ext"
	# inparalog_int (from the same species and the same orthology group)
	evtype[ ( eva["ev_type"] == "D" ) & ( eva["cluster_name_x"] == eva["cluster_name_y"] ) & ( eva["in_sps"] == eva["out_sps"] ) ] = "inparalog_int"
	# inparalog_ext (from the same species and a different orthology group)
	evtype[ ( eva["ev_type"] == "D" ) & ( eva["cluster_name_x"] != eva["cluster_name_y"] ) & ( eva["in_sps"] == eva["out_sps"] ) ] = "inparalog_ext"
	eva["ev_type"] = evtype

	return eva


#####################
####### MAIN ########
#####################

if __name__ == '__main__':

	# read phylogeny, find speciation events, create network, do clustering
	evs, eva, phy, phy_lis, clu = parse_phylo(phy_fn=phy_fn, phy_id=phy_id, do_allpairs=do_allpairs, do_root=do_root, clusters_function_string=clusters_function_string, do_sps_reconciliation=do_sps_reconciliation)
	
	# make human readable cluster names (instead of integers)
	clu["cluster_name"] = ogprefix + clu["cluster"].astype(str)
	clu["node_ref"] = ""
	clu["node_ref_support"] = ""
	
	# find cluster supports (support in oldest node in cluster)
	clu["support"] = find_support_cluster(clu=clu, phy=phy, cluster_label="cluster")
	
	# infer gene-to-gene orthology relationship classes (inparalog/outparalog/ortholog)
	if do_allpairs:
		
		eva = annotate_event_type(eva=eva, clu=clu, clutag="cluster_name")
		eva = eva[["in_gene","out_gene","ev_type"]]
		eva = pd.DataFrame(eva).dropna()
	
	# try to annotate reference sequences
	if do_ref:
		
		# load ref
		ref = pd.read_csv(ref_fn, sep="\t", names=["gene","name"])
		
		# remove ref nodes not in phylogeny
		ref =  ref[ np.isin(element=ref["gene"].values, test_elements=phy_lis) ]
		
		# report which reference sequences can be found within cluster
		clu["cluster_ref"]  = ref_tagcluster(clu=clu, evs=evs, ref=ref, ref_spi=refsps, label_if_no_annot="NA", clean_gene_names=clean_gene_names)
		clu["cluster_name"] = ogprefix + clu["cluster"].astype(str) + ":" + clu["cluster_ref"].astype(str)
		print_attributes    = ["cluster_name"]
		
		# find named orthologs anywhere in the phylogeny
		clu["node_ref"], clu["node_ref_support"] = ref_known_any(clu=clu, evs=evs, ref=ref, ref_spi=refsps, syn_nod=None)
		print_attributes.append("node_ref")
		
		# extend cluster-wise annotations
		if min_support_transfer is not None:
			
			# extend cluster-specific labels
			clu["extended_clusters"], clu["extended_labels"] = find_close_monophyletic_clusters(clu=clu, phy=phy, ref_label="cluster_ref", ref_NA_label="NA", cluster_label="cluster", min_support_transfer=min_support_transfer)
			ixs_to_rename = np.where(clu["extended_clusters"].values != None)[0]
			clu.loc[ ixs_to_rename, "cluster_name" ] = ogprefix + clu.loc[ ixs_to_rename, "cluster" ].astype(str) + ":like:" + clu.loc[ ixs_to_rename, "extended_labels" ].astype(str)
			
	else:
		
		print_attributes = ["cluster_name"]		
	
	# print phylogeny
	write_tree(phy=phy, out="%s/%s.ortholog_groups.newick" % (out_fn,phy_id), evc=clu, attributes=print_attributes, sep=" | ", do_print=do_print, cut_gene_names=cut_gene_names)
	
	# save clusters
	clu_print = clu[["node","cluster_name","support","node_ref","node_ref_support"]]
	clu_print.columns = ["gene","orthogroup","orthogroup_support","reference_ortholog","reference_support"]
	clu_print.to_csv("%s/%s.ortholog_groups.csv" % (out_fn,phy_id), sep="\t", index=None, mode="w")
	
	# save gene pair information (orthologs and all genes)
	evs.to_csv("%s/%s.pairs_orthologs.csv" %  (out_fn,phy_id), sep="\t", index=None, mode="w")
	if do_allpairs:
		eva.to_csv("%s/%s.pairs_all.csv" %  (out_fn,phy_id), sep="\t", index=None, mode="w")
	
	logging.info("%s Done" % phy_id)

