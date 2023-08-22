import ete3
import numpy as np
import pandas as pd
import argparse

# argument parser
arp = argparse.ArgumentParser()

# Add the arguments to the parser
arp.add_argument("-tree",  "--tree",  required=True,                        help="newick file with species tree")
arp.add_argument("-ort",   "--ort",   required=True,                        help="table of orthologs (long format, with headers; min info: OG <tab> gene1)")
arp.add_argument("-out",   "--out",   required=True,                        help="output prefix")
arp.add_argument("-gcol",  "--gcol",  required=False, default="gene",       help="OPTIONAL: name of gene column in --ort/-t table. Default is \"gene\"")
arp.add_argument("-ccol",  "--ccol",  required=False, default="orthogroup", help="OPTIONAL: name of orthogroup column in --ort/-t table. Default is \"orthogroup\"")
arp.add_argument("-split", "--split", required=False, default="_",          help="OPTIONAL: character to split species and sequence names. Default is \"_\", e.g. Human_genename. WARNING: use quotation marks, e.g. -split \"_\" or -split \"|\"")
arl = vars(arp.parse_args())

# input variables
phs_fn   = arl["tree"]
ort_fn   = arl["ort"]
out_fn   = arl["out"]
gene_col = arl["gcol"]
clus_col = arl["ccol"]
split_ch = arl["split"].replace("\"","")

# out_fn = "ara"
# phs_fn = "species.newick"
# ort_fn = "tables_grouped/phy.HLH.ortholog_groups.csv"
# ref_sps = "Hsap"
# gene_col = "gene"
# clus_col = "orthogroup"
# split_ch = "_"

# find nodes at which each species pair is gained
def do_species_orig_dict(phs):

	sps_list = phs.get_leaf_names()
	nod_list = [ i for i in phs.get_descendants() ]
	anc_list = [ phs.name ] + [ i.name for i in nod_list if not i.is_leaf() ]
	anc_list = anc_list[::-1]

	species_orig_dict = dict()
	species_ages_dict = dict()

	for n,i in enumerate(sps_list):
		species_orig_dict[i] = dict()
		species_ages_dict[i] = dict()
		for m,j in enumerate(sps_list):
			if n != m:

				# find common ancestor
				e = phs.get_common_ancestor(i,j)
				
				# find its name (or create one from descendant nodes if no name given)
				if e.name == "" :
					species_orig_dict[i][j] = e.get_leaf_names(topology_only = True)[0] + "," + e.get_leaf_names(topology_only = True)[-1]
				else:
					species_orig_dict[i][j] = e.name

				# find its age
				species_ages_dict[i][j] = e.get_farthest_leaf(topology_only = True)[1]

			elif n == m:
				species_orig_dict[i][j] = i
				species_ages_dict[i][j] = 0

	return species_orig_dict, species_ages_dict, sps_list, anc_list

def collapsed_leaf(node):
    if len(node2labels[node]) == 1:
       return True
    else:
       return False

# loop through list of orthogroups and calculate ages
def do_ancestral_reconstruction(ort, phs, species_orig_dict, species_ages_dict, sps_list, anc_list, clus_col=clus_col, gene_col=gene_col, split_ch=split_ch):

	# list of orthoclusters
	clus_lis = np.unique(ort[clus_col])
	
	# empty lists for outputs
	gain_lis = np.zeros(len(clus_lis), dtype=object)
	loss_lis = np.zeros(len(clus_lis), dtype=object)
	pres_lis = np.zeros(len(clus_lis), dtype=object)
	
	# list of nodes (extant and ancestral)
	cols_lis = sps_list + anc_list

	dat = pd.DataFrame()
	mat_gain = pd.DataFrame(0 , columns = cols_lis , index=clus_lis)
	mat_loss = pd.DataFrame(0 , columns = cols_lis , index=clus_lis)
	mat_pres = pd.DataFrame(0 , columns = cols_lis , index=clus_lis)

	# loop for each cluster in dataset
	for n,c in enumerate(clus_lis):

		# genes present in cluster
		nod_clu = ort[ort[clus_col] == c][gene_col].values            
		# species present in cluster  
		sps_clu = np.unique([ m.split(split_ch)[0] for m in nod_clu ])
		sps_clu_string = ",".join(sorted(sps_clu))
		# ... and their counts
		sps_clu_c = np.unique([ m.split(split_ch)[0] for m in nod_clu ], return_counts=True)[1]


		### GAINS ###
		# find origin of cluster
		# if present in only one species, it's a singleton
		# if present in more than one species, go fetch their last common ancestor
		if len(sps_clu) == 1:
			clu_gain = sps_clu[0]
		else:
			age = -1 # start with negative age to accommodate species-specific (age=0)
			for ni,si in enumerate(sps_clu):
				for nj,sj in enumerate(sps_clu):
					if ni < nj:
						age_ij = species_ages_dict[si][sj] # divergence time of current species pair
						if age_ij > age: # if divergence time of current pair is older than age, reassign age
							age = int(age_ij)
							clu_gain = species_orig_dict[si][sj]

		# store cluster gain matrix
		mat_gain[clu_gain][c] = 1
		# store per-og matrix
		gain_lis[n] = clu_gain
		pres_lis[n] = sps_clu_string
		
		### PRESENCE ###
		# a cluster is present in a node of the species tree if it is present in
		# any of the intermediate nodes between extant nodes and the node of origin
		phs_gain_node     = phs.search_nodes(name=clu_gain)[0]
		phs_gain_node_sps = phs_gain_node.get_leaf_names()         # find all descendant species
		set_pres = set()
		for ni,si in enumerate(sps_clu):
			phs_si = phs.search_nodes(name=si)[0]
			while True:
				set_pres.add(phs_si.name)
				if phs_si == phs_gain_node:
					break
				else:
					phs_si = phs_si.up

		# fill matrix with presence data
		for clu_pres in set_pres:
			mat_pres[clu_pres][c] = 1
		# modify extant species to contain counts rather than 1/0
		for clu_prei,clu_pres in enumerate(sps_clu):
			mat_pres[clu_pres][c] = sps_clu_c[clu_prei]


		### LOSSES ###
		# cluster is lost in a node of the species tree if it is absent in this node and in all its descendants 
		# init a set of losses with all species where the cluster is missing (starting from the descendants of the origin node)
		set_loss_init = {element for element in set(phs_gain_node_sps) if element not in set(sps_clu)}
		set_loss = set()
		for ni,si in enumerate(set_loss_init):
			phs_si = phs.search_nodes(name=si)[0]
			while True:

				# get species hanging from this node
				phs_si_sps = phs_si.get_leaf_names()

				# if the descendant species are a subset of the species that lack this cluster,
				# carry on (go further up in the tree), we have not reached the loss node yet.
				# if any of the descendant species lacks the cluster, add the the node from the 
				# previous iteration to the the set of loss nodes and break the loop
				if set(phs_si_sps).issubset(set_loss_init):
					phs_candidate = phs_si
					phs_si = phs_si.up
				else:
					set_loss.add(phs_candidate.name)
					break

		# store losses in matrix
		for clu_loss in set_loss:
			mat_loss[clu_loss][c] = 1
		# store per-og matrix
		loss_lis[n] = ",".join(sorted(set_loss))

		### LOG  ###
		print("# %i/%i | %s | %s | %s" % (n+1,len(clus_lis),c, clu_gain,sps_clu_string))
		# print("present:", set_pres)
		# print("losses: ", set_loss)


	# prepare output
	dat = pd.DataFrame({
		"orthogroup" :  clus_lis,
		"presence":     pres_lis,
		"gain" :        gain_lis,
		"loss" :        loss_lis,
		"n_gains" :     mat_gain.sum(axis=1),
		"n_losses" :    mat_loss.sum(axis=1),
		"n_presences" : mat_pres.sum(axis=1) 
	}, columns=[ "orthogroup","presence","gain","loss","n_gains","n_losses","n_presences"] )

	# output
	return dat, mat_gain, mat_loss, mat_pres


#### MAIN WORK ####

# load input tree
print("# Load tree from %s" % phs_fn)
phs = ete3.PhyloTree("%s" % (phs_fn), format=1)
# assign species names to tree
phs.set_species_naming_function(lambda node: node.name)
# resolve polytomies in a random fashion
#phs.resolve_polytomy(recursive=True)

# load orthoclusters
print("# Load orthoclusters from %s" % ort_fn)
ort = pd.read_csv(ort_fn, sep="\t")
ort = ort[[gene_col,clus_col]]

# obtain species-to-species dictionary of relative ages
print("# Species-to-species relative ages from %s" % phs_fn)
species_orig_dict,species_ages_dict, sps_list, anc_list = do_species_orig_dict(phs=phs)
dat,mat_gain,mat_loss,mat_pres = do_ancestral_reconstruction(ort=ort,phs=phs, species_orig_dict=species_orig_dict, species_ages_dict=species_ages_dict, sps_list=sps_list, anc_list=anc_list)

# save output: per orthogroup
print("# Save output")
dat.to_csv("%s.dollo_summary.csv" % out_fn, sep="\t", index=None, mode="w")
mat_gain.to_csv("%s.dollo_gain.csv" % out_fn, sep="\t",  mode="w")
mat_loss.to_csv("%s.dollo_loss.csv" % out_fn, sep="\t", mode="w")
mat_pres.to_csv("%s.dollo_pres.csv" % out_fn, sep="\t", mode="w")

print("# Fet!")


