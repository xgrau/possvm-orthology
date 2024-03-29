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
def parse_phylo(phy_fn, phy_id, do_root, do_allpairs, outgroup=outgroup):

	# load input
	phy = ete3.PhyloTree("%s" % (phy_fn))
	logging.info("%s num nodes = %i" % (phy_id,len(phy)))

	# assign species names to tree
	phy.set_species_naming_function(lambda node: node.name.split(split_ch)[0] )

	# resolve polytomies (randomly)
	phy.resolve_polytomy(recursive=True)
	
	# try to find root if unrooted
	if do_root:
	
		# shall we do it with iterative midpoint rooting?
		if itermidroot is not None:

			# niter = max(int(len(phy)/10),10)
			# niter = int(len(phy)*0.9)
			niter = itermidroot
			num_evs_per_iter = np.zeros(niter)
			out_nod_per_iter = np.empty(niter,	dtype=object)

			# then, iterate to try to find better candidates
			phy_it = phy
			phy_outgroup_it = phy_it.get_midpoint_outgroup()
			for roi in range(niter):

				# in subsequent iterations, assign zero distance to the previous outgroup edge, and try to find second-best candidate
				phy_outgroup_it.dist = 0.0
				phy_outgroup_it = phy_it.get_midpoint_outgroup()
				phy_it.set_outgroup(phy_outgroup_it)

				# parse events and re-do clustering
				evs_it, _, _, phy_lis_it = parse_events(phy=phy_it, outgroup=outgroup, do_allpairs=False, min_support_node=min_support_node)
				clu_it = clusters_mcl(evs=evs_it, node_list=phy_lis_it, verbose=False)

				# store number of orthogroups in this particular iteration
				num_evs_per_iter[roi] = len(np.unique(clu_it["cluster"].values))
				out_nod_per_iter[roi] = phy_outgroup_it

				logging.info("%s Iterative midpoint root | %i/%i | n OGs = %i" % (phy_id,roi+1,niter,num_evs_per_iter[roi]))

			# select outgroup that minimises number of OGs (more agglomerative)
			phy_outgroup_ix = np.argmin(num_evs_per_iter)
			phy_outgroup = out_nod_per_iter[phy_outgroup_ix]

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
	evs, eva, phy, phy_lis = parse_events(phy=phy, outgroup=outgroup, do_allpairs=do_allpairs, min_support_node=min_support_node)
	clu = clusters_mcl(evs=evs, node_list=phy_lis)

	# output from event parsing
	return evs, eva, phy, phy_lis, clu


# parse phylogenies with ETE to obtain a network-like table defining 
# orthologous relationships, using the species overlap algorithm
def parse_events(phy, outgroup, do_allpairs, min_support_node=0):

	# list of genes in phylogeny
	phy_lis = phy.get_leaf_names()

	# find evolutionary events (duplications and speciations)
	evev = phy.get_descendant_evol_events(sos_thr=sos)

	# speciation events
	evs    = np.empty((len(evev)*len(evev), 5), dtype="object")
	evs[:] = np.nan
	n = 0
	for ev in evev:
		# if ev.etype == "S" :   # normally, this would work. But to avoid splitting monophyletic single species clades into different clusters, we'll add single-species monophyletic paralogs as orthologs
		# retrieve list of sps in ingroup and outgroup:
		sps_in = np.unique([ i.split(split_ch)[0] for i in ev.in_seqs ])
		sps_ou = np.unique([ i.split(split_ch)[0] for i in ev.out_seqs ])
		# check if node is a speciation node, or a duplication node where both descendant branches have exactly one species, and this is the same species
		if (ev.etype == "S" or (ev.etype == "D" and len(sps_in) == len(sps_ou) == 1 and sps_in == sps_ou)) and ev.branch_supports[0] >= min_support_node:
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

	# drop outgroup species, if any
	if len(outgroup) > 0:
		in_evs_sps =  [ i.split(split_ch)[0] in set(outgroup) for i in evs["in_gene"]  ]
		out_evs_sps = [ i.split(split_ch)[0] in set(outgroup) for i in evs["out_gene"] ]
		evs = evs.drop(np.where(np.logical_or(in_evs_sps, out_evs_sps))[0])
		phy_lis =  [ i for i in phy_lis if i.split(split_ch)[0] not in set(outgroup) ]

	if do_allpairs:
		# duplications and speciation events
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


# function to cluster a network-like table of orthologs (from ETE) 
def clusters_lpa(evs, node_list, cluster_label="cluster", label_if_no_annot="NA"):

	# clustering: create network
	logging.info("Create network")
	evs_e = evs[["in_gene","out_gene","branch_support"]]
	evs_n = nx.convert_matrix.from_pandas_edgelist(evs_e, source="in_gene", target="out_gene", edge_attr="branch_support")
	evs_n.add_nodes_from(node_list)
	
	# clustering: asynchronous label propagation
	logging.info("Find communities LPA")
	clu_c = community.asyn_lpa_communities(evs_n, seed=11)and min_support_node > min_support_node
	clu_c = { frozenset(c) for c in clu_c }
	logging.info("Find communities LPA num clusters = %i" % len(clu_c))
	clu_c_clu = [ i for i, cluster in enumerate(clu_c) for node in cluster ]
	clu_c_noi = [ node for i, cluster in enumerate(clu_c) for node in cluster ]

	# clustering: save output
	clu = pd.DataFrame( { 
		"node"    :  clu_c_noi,
		cluster_label : clu_c_clu,
	}, columns=["node",cluster_label])
	logging.info("Find communities LPA | num clustered genes = %i" % len(clu))

	return clu

# function to cluster a network-like table of orthologs (from ETE) 
def clusters_louvain(evs, node_list, cluster_label="cluster", label_if_no_annot="NA"):

	# clustering: create network
	logging.info("Create network")
	evs_e = evs[["in_gene","out_gene","branch_support"]]
	evs_n = nx.convert_matrix.from_pandas_edgelist(evs_e, source="in_gene", target="out_gene", edge_attr="branch_support")
	evs_n.add_nodes_from(node_list)

	# clustering: Louvain
	import community as community_louvain

	logging.info("Find communities Louvain")
	clu_x = community_louvain.best_partition(evs_n)
	clu_c = {}
	for k, v in clu_x.items():
		clu_c[v] = clu_c.get(v, [])
		clu_c[v].append(k)
	logging.info("Find communities Louvain | num clusters = %i" % len(clu_c))
	clu_c_noi = [ n for i,c in enumerate(clu_c) for n in clu_c[c] ]
	clu_c_clu = [ c for i,c in enumerate(clu_c) for n in clu_c[c] ]

	# clustering: save output
	clu = pd.DataFrame( { 
		"node"    :  clu_c_noi,
		cluster_label : clu_c_clu,
	}, columns=["node",cluster_label])
	logging.info("Find communities Louvain | num clustered genes = %i" % len(clu))

	return clu


# function to cluster a network-like table of orthologs (from ETE) 
# using MCL
def clusters_mcl(evs, node_list, inf=1.5, verbose=True):

	if len(evs) > 0:

		# MCL clustering: create network
		if verbose:
			logging.info("Create network")
		evs_e = evs[["in_gene","out_gene","branch_support"]]
		evs_n = nx.convert_matrix.from_pandas_edgelist(evs_e, source="in_gene", target="out_gene", edge_attr="branch_support")
		evs_n.add_nodes_from(node_list)

		evs_n_nodelist = [ node for i, node in enumerate(evs_n.nodes()) ]
		evs_m = nx.to_scipy_sparse_matrix(evs_n, nodelist=evs_n_nodelist)
		# MCL clustering: run clustering
		if verbose:
			logging.info("MCL clustering, inflation = %.3f" % (inf))
		mcl_m  = markov_clustering.run_mcl(evs_m, inflation=inf, pruning_threshold=0) #why pruning threshold HAS to be zero?
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

		mcl_c_clu = [ i for i in range(len(node_list)) ]

	# output
	clu = pd.DataFrame( { 
		"node"    : mcl_c_nod,
		"cluster" : mcl_c_clu,
	}, columns=["node","cluster"])
	if verbose:
		logging.info("MCL clustering, num clustered genes = %i" % (len(clu)))

	return clu


# function to cluster a network-like table of orthologs (from ETE) 
def clusters_dist(phy, k, cluster_label="cluster"):

	logging.info("Obtain pairwise distances")

	# list of genes in phylogeny
	phy_lis = phy.get_leaf_names()

	# distance matrix
	dis_m = np.zeros(shape=((len(phy_lis), len(phy_lis))))
	for n,i in enumerate(phy_lis):
		for m,j in enumerate(phy_lis):
			if n < m:
				dis_m[n,m] = phy.get_distance(i,j)

	# only upper triangle is populated, but we can transpose and sum to get lower triangle too
	dis_m = dis_m + dis_m.T
	dis_d = 1/dis_m

	# graph from matrix
	dis_g = nx.convert_matrix.from_numpy_matrix(dis_d)
	# clu_c = community.asyn_lpa_communities(dis_g, weight="weight")
	clu_c = community.asyn_fluidc(dis_g, k=k)

	# # create reducer
	# reducer = umap.UMAP()
	# embedding = reducer.fit_transform(dis_m)

	# plt.scatter( embedding[:, 0], embedding[:, 1] , c=clu["cluster"].values) 
	# plt.show()

	clu_c = {frozenset(c) for c in clu_c}
	logging.info("clustering, num clusters = %i" % len(clu_c))
	clu_c_clu = [ i for i, cluster in enumerate(clu_c) for node in cluster ]

	return clu_c_clu



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

def ref_annot(clu, evs, ref_spi, syn_nod, label_if_no_annot=""):

	logging.info("Add annotations: orthology to genes from reference species")

	evs["in_sps"]  = evs["in_gene"].apply(lambda c: c.split(split_ch)[0])
	evs["out_sps"] = evs["out_gene"].apply(lambda c: c.split(split_ch)[0])

	cluster_ref = []
	for noi in clu["node"]:

		# which cluster?
		c = clu[clu["node"] == noi]["cluster"].values[0]

		# find reference sequences WITHIN THIS cluster 
		# and create a dictionary to alphabetic short codes
		ref_nodes = clu[(clu["sps"] == ref_spi) & (clu["cluster"] == c)]["node"].values
		ref_codes = list(string.ascii_letters[0:len(ref_nodes)])
		ref_dicti = dict()
		for m,r in enumerate(ref_nodes):
			ref_dicti[r] = ref_codes[m]

		# find if gene is orthologous to any ref sequences
		r1 = evs[(evs["in_gene"] == noi) & (evs["out_sps"] == ref_spi)]["out_gene"].values
		r2 = evs[(evs["out_gene"] == noi) & (evs["in_sps"] == ref_spi)]["in_gene"].values
		ra = np.unique(np.concatenate((r1,r2)))
		ra_is_in_cluster = np.isin(element=ra, test_elements=ref_nodes)
		ra = ra[ra_is_in_cluster]

		# if reference sequences in ra have synonyms and these synonyms are refs too, add synonyms to ra
		# ra = add_synonymous_nodes(nodes=ra, syn_nod=syn_nod, ref_nodes=ref_nodes)

		# use cluster-specific alphabetic codes for ref sequences
		rc = np.unique(np.sort([ ref_dicti[r] for r in ra ]))
		rc = ''.join(rc)

		# name of sps-specific cluster: cluster_0ab_sps
		cluster_ref.append(str(c)+rc)

	ix_ref_spi      = np.where((clu["sps"] == ref_spi).values)
	cluster_ref_table = pd.DataFrame( { 
		"node" :        clu["node"].values[ix_ref_spi],
		"cluster" :     clu["cluster"].values[ix_ref_spi],
		"cluster_ref" : np.array(cluster_ref)[ix_ref_spi]
	} , columns=["node","cluster","cluster_ref"])

	return cluster_ref, cluster_ref_table

# if any node in an array has synonyms, add synonyms to array
def add_synonymous_nodes(nodes, syn_nod, ref_nodes=None):
	
	if ref_nodes is None:
		ref_nodes = nodes
	
	for r in nodes:
		for syn in syn_nod:
			for s in syn:
				if r in syn and np.isin(s, ref_nodes) and s is not r:
					nodes = np.append(nodes, s)
	return nodes

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

		# if reference sequences in ra have synonyms and these synonyms are refs too, add synonyms to ra
		# ra = add_synonymous_nodes(nodes=ra, syn_nod=syn_nod, ref_nodes=ref_nodes)

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


def find_close_clusters(clu, phy, extension_ratio_threshold, ref_label="cluster_ref", ref_NA_label="NA", cluster_label="cluster"):
	
	logging.info("Add annotations: extend annotations to close groups | extension ratio threshold = %.3f" % extension_ratio_threshold)

	# input and output lists
	has_label_ix = np.where(clu[ref_label] != ref_NA_label)[0]
	has_no_label_ix = np.where(clu[ref_label] == ref_NA_label)[0]
	clusters_labeled = np.unique(clu[cluster_label][has_label_ix].values)
	clusters_nonlabeled = np.unique(clu[cluster_label][has_no_label_ix].values)
	
	# extended annotations: init arrays with previous information
	extended_annots = clu[ref_label].values
	extended_clusters = clu[cluster_label].values
	# extended_annots = np.repeat("NA", clu[ref_label].shape[0])
	# extended_clusters = np.repeat(0, clu[ref_label].shape[0])
	num_extended = 0

	if len(clusters_labeled) > 1:

		# loop for each non-labeled cluster
		for ci in clusters_nonlabeled:

			# list of nodes in non-labeled cluster ci
			nodes_i = clu[clu[cluster_label] == ci]["node"].values.tolist()

			# index vector
			nodes_in_cluster_ci_ix = np.where(clu[cluster_label] == ci)[0]
			
			# if there is only one node, root is the same node (otherwise defaults to tree root!)
			if len(nodes_i) > 1:
				root_i  = phy.get_common_ancestor(nodes_i)
			else:
				root_i = nodes_i[0]

			# init values for cloest node search
			first_closest_d_phyl = 1e6
			first_closest_c = ci
			second_closest_d_phyl = 1e6
			second_closest_c = ci

			# loop labeled clusters to find closest
			for cj in clusters_labeled:

				nodes_j = clu[clu[cluster_label] == cj]["node"].values.tolist()
				if len(nodes_j) > 1:
					root_j = phy.get_common_ancestor(nodes_j)
				else:
					root_j = nodes_j[0]

				# calculate topological and phylogenetic distances between roots of clusters ci and cj
				dist_ij_phyl = phy.get_distance(root_i, root_j)

				# check if current labeled root (cj) is closest to unlabeled root (ci)
				if dist_ij_phyl < first_closest_d_phyl:
					first_closest_d_phyl = dist_ij_phyl
					first_closest_c = cj

			# same for second closest
			for cj in clusters_labeled:

				nodes_j = clu[clu[cluster_label] == cj]["node"].values.tolist()
				if len(nodes_j) > 1:
					root_j = phy.get_common_ancestor(nodes_j)
				else:
					root_j = nodes_j[0]

				# calculate topological and phylogenetic distances between roots of clusters ci and cj
				dist_ij_phyl = phy.get_distance(root_i, root_j)

				# check if current labeled root (cj) is second closest to unlabeled root (ci)
				if dist_ij_phyl < second_closest_d_phyl and dist_ij_phyl > first_closest_d_phyl:
					second_closest_d_phyl = dist_ij_phyl
					second_closest_c = cj

			# add min value to avoid zero divisions
			first_closest_d_phyl = np.max((first_closest_d_phyl, 1e-6))

			# distance ratio (second to first)
			d2d1 = second_closest_d_phyl / first_closest_d_phyl

			if d2d1 > extension_ratio_threshold:

				# which is the annotation in the closest group?
				first_closest_annot = np.unique(
					clu[ref_label] [ (clu[cluster_label] == first_closest_c) & (clu[ref_label] != ref_NA_label) ].values
				)

				# assign new label
				extended_annots[nodes_in_cluster_ci_ix] = first_closest_annot
				extended_clusters[nodes_in_cluster_ci_ix] = first_closest_c

				# counter
				num_extended = num_extended + 1

				# print("# OG%i will now be extOG%i || d2/d1=%.3f || d1 OG%i = %.3f || d2 OG%i = %.3f" % (ci, first_closest_c, d2d1, first_closest_c, first_closest_d_phyl, second_closest_c, second_closest_d_phyl))

			else:
				
				pass
				# print("# OG%i will remain as is || d2/d1=%.3f || d1 OG%i = %.3f || d2 OG%i = %.3f" % (ci, d2d1, first_closest_c, first_closest_d_phyl, second_closest_c, second_closest_d_phyl))



	logging.info("Add annotations: extend annotations to close groups | %i labels transferred" % num_extended)
	
	return extended_clusters, extended_annots


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


# find clades of monophyletic sequences (useful for seq annotation)
def find_monophyletic_expansion_refsps(phy, ref_sps):

	syn_sets = list()
	for ref_spi in ref_sps:

		logging.info("Find synonyms: monophyletic | sps = %s" % (ref_spi))

		# list of species and gene names of descending leaves of each node in the tree
		node_species = phy.get_cached_content(store_attr="species")
		node_geneids = phy.get_cached_content(store_attr="name")

		# find sets of monophyletic sequences from the ref sps
		for n,i in enumerate(node_species):
			if (len(node_species[i]) == 1) and (ref_spi in node_species[i]) and (len(list(node_geneids[i])) > 1):
				syn_sets.append(node_geneids[i])

	# remove redundancy in these sets (check for overlapping sets and merge into larger sets)
	syn_sets_nr = syn_sets
	for n,_ in enumerate(syn_sets_nr):
		for m,_ in enumerate(syn_sets_nr):
			if n != m:
				if len(syn_sets_nr[n].intersection(syn_sets_nr[m])) > 1:
					syn_sets_nr[n] =  syn_sets_nr[n].union(syn_sets_nr[m]) 
	syn_sets_nr = np.unique(syn_sets_nr)

	logging.info("Find synonyms: monophyletic | n = %i syn sets found" % len(syn_sets_nr))

	return syn_sets_nr

# find groups of genes from the same species that are very similar
# according to phyloenetic distance, and output set of synonyms
def find_close_gene_pairs_refsps(phy, ref_sps, dist_thr=0.0, split_ch=split_ch):

	syn_sets = list()
	for ref_spi in ref_sps:

		logging.info("Find synonyms: close paralogs | sps = %s | d = %.2f" % (ref_spi, dist_thr))

		list_genes = phy.get_leaf_names()
		list_species = [ c.split(split_ch)[0] for c in list_genes ]

		for i in phy.iter_leaves():
			if i.species == ref_spi:
				for m,j in enumerate(list_genes):
					if i.name != j and i.species == list_species[m] and i.get_distance(j) < dist_thr:
						syn_sets.append(set((i.name, j)))
	
	# remove redundancy in these sets (check for overlapping sets and merge into larger sets)
	syn_sets_nr = syn_sets
	for n,_ in enumerate(syn_sets_nr):
		for m,_ in enumerate(syn_sets_nr):
			if n != m:
				if len(syn_sets_nr[n].intersection(syn_sets_nr[m])) > 1:
					syn_sets_nr[n] = syn_sets_nr[n].union(syn_sets_nr[m])

	syn_sets_nr = np.unique(syn_sets_nr)

	logging.info("Find synonyms: close paralogs | n = %i syn sets found" % len(syn_sets_nr))

	return syn_sets_nr


def sanitise_genename_string(string, splitstring="/"):

	import re
	
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

