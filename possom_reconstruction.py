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
					species_orig_dict[i][j] = e.get_leaf_names()[0] + "," + e.get_leaf_names()[-1]
				else:
					species_orig_dict[i][j] = e.name

				# find its age
				species_ages_dict[i][j] = e.get_farthest_leaf()[1]

			elif n == m:
				species_orig_dict[i][j] = i
				species_ages_dict[i][j] = 0

	return species_orig_dict, species_ages_dict, sps_list, anc_list


# loop through list of orthogroups and calculate ages
def do_ancestral_reconstruction(ort, phs, species_orig_dict, species_ages_dict, sps_list, anc_list, clus_col=clus_col, gene_col=gene_col, split_ch=split_ch):

	# list of orthoclusters
	clus_lis = np.unique(ort[clus_col])
	
	# empty lists for outputs
	gain_lis = np.zeros(len(clus_lis), dtype=object)
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
		# ... AND THEIR COUNTS
		sps_clu_c = np.unique([ m.split(split_ch)[0] for m in nod_clu ], return_counts=True)[1]


		### GAINS ###
		# find origin of cluster
		a = -1
		for ni,si in enumerate(sps_clu):
			for nj,sj in enumerate(sps_clu):
				if ni<=nj:
					t = species_ages_dict[si][sj] # divergence time of current species pair
					if t > a: # if divergence time of current pair is older than age, reassign age (a) and named ages (o)
						a = int(t)
						clu_gain = species_orig_dict[si][sj]

		# store cluster gain matrix
		mat_gain[clu_gain][c] = 1
		# store per-og matrix
		gain_lis[n] = clu_gain
		pres_lis[n] = sps_clu_string
		
		# log a bit
		print("# %i / %i | %s | %s | %s" % (n+1,len(clus_lis),c, clu_gain,sps_clu_string))

		### LOSSES ###
		# find which species should have this cluster but don't
		phs_o = phs.search_nodes(name=clu_gain)[0] # find origin node in species tree
		phs_o_sps = phs_o.get_leaf_names()         # find all descendant species
		sps_miss = list({element for element in set(phs_o_sps) if element not in set(sps_clu)}) # missing species (from origin node only)

		# if cluster is missing from >=1 species, try to agglomerate losses to shared ancestral nodes
		if len(sps_miss) > 1:

			sps_miss_iter = set(sps_miss[:])
			sps_miss_init = set(sps_miss[:])
			
			while True:

				sps_miss_loss = set()
				for ni,si in enumerate(sps_miss_iter):

					# find sister node
					phs_si = phs.search_nodes(name=si)[0]
					phs_si_sister = phs_si.get_sisters()[0]
					loss_in_node = si # init loss to same species

					# look at all other missing species and check if we can take loss back to a shared ancestor
					for nj,sj in enumerate(sps_miss_iter):
						if ni != nj:
							phs_sj = phs.search_nodes(name=sj)[0]
							# if sister node is equals another missing species, move loss to shared ancestor node:
							if phs_si_sister == phs_sj:
								loss_in_node = phs.get_common_ancestor(phs_si,phs_sj).name

					# store loss nodes
					sps_miss_loss.add(loss_in_node)

				if sps_miss_loss == sps_miss_init:
					break
				else:
					sps_miss_iter = sps_miss_loss
					sps_miss_init = sps_miss_loss
		
		# if cluster is missing in less than one species, two options:
		# 1. it either got lost in that species
		# 2. it's never been lost
		# ergo, we'll use an empty losses set
		else:
			sps_miss_loss = set(sps_miss)

		# store losses in matrix
		for clu_loss in sps_miss_loss:
			mat_loss[clu_loss][c] = 1


		### PRESENCE ###
		# first, fill matrix with extant presences
		# use extant counts instead of simply "1"
		for clu_prei,clu_pres in enumerate(sps_clu):
			mat_pres[clu_pres][c] = sps_clu_c[clu_prei]

		# second-i, find ancestral presence nodes
		# ancestral presence = descendant nodes from the original gain that are not losses themselves
		phs_gain_node = phs.search_nodes(name=clu_gain)[0]
		descendants_from_gain_node = phs_gain_node.get_descendants()
		descendants_from_gain_node = [ i.name for i in descendants_from_gain_node if not i.is_leaf() ]
		presents_from_gain_node = list({element for element in set(descendants_from_gain_node) if element not in set(sps_miss_loss)})

		# second-ii, fill matrix with ancestral presences
		for clu_pres in presents_from_gain_node:
			mat_pres[clu_pres][c] = 1


	# prepare output
	dat = pd.DataFrame({
		"orthogroup" :  clus_lis,
		"presence":     pres_lis,
		"gain" :        gain_lis,
		"n_gains" :     mat_gain.sum(axis=1),
		"n_losses" :    mat_loss.sum(axis=1),
		"n_presences" : mat_pres.sum(axis=1) 
	}, columns=[ "orthogroup","presence","gain","n_gains","n_losses","n_presences"] )

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


