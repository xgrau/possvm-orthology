import ete3
import numpy as np
import pandas as pd
import argparse

# argument parser
arp = argparse.ArgumentParser()

# Add the arguments to the parser
arp.add_argument("-t", "--tree",  required=True,                        help="newick file with species tree")
arp.add_argument("-i", "--in",    required=True,                        help="table of orthologs (long format, with headers; min info: OG <tab> gene1)")
arp.add_argument("-o", "--out",   required=True,                        help="output prefix")
arp.add_argument("-g", "--gcol",  required=False, default="gene",       help="OPTIONAL: name of gene column in --in/-i table. Default is \"node\"")
arp.add_argument("-c", "--ccol",  required=False, default="orthogroup", help="OPTIONAL: name of orthogroup column in --in/-i table. Default is \"og_cluster\"")
arp.add_argument("-s", "--split", required=False, default="_",          help="OPTIONAL: character to split species and sequence names. Default is \"_\", e.g. Human_genename. WARNING: use quotation marks, e.g. -split \"_\" or -split \"|\"")
arl = vars(arp.parse_args())

# input variables
phs_fn   = arl["tree"]
ort_fn   = arl["in"]
out_fn   = arl["out"]
gene_col = arl["gcol"]
clus_col = arl["ccol"]
split_ch = arl["split"].replace("\"","")


# define a dictionary of species-two-species relative ages
# for each species in the phyolgeny
def species_age_dict(phs):
		
	# init dict of dicts
	sps_age_dict = dict()
	sps_out_dict = dict()
	sps_nod_dict = dict()
	sps_list     = phs.get_leaf_names()
	age_root     = 0
	for n,i in enumerate(sps_list):

		# init dict of ages for species i
		sps_age_dict[i] = dict()
		for m,j in enumerate(sps_list):
			if n != m:
				e = phs.get_common_ancestor(i,j) # define species subtree between species i and j
				a = e.get_farthest_leaf()[1]     # distance to the farthest leaf in subtree e (age)
			if n == m:
				a = 0

			# new entry j in dict for species i
			sps_age_dict[i][j] = int(a)

			# is this the most ancient node?
			if sps_age_dict[i][j] > age_root:
				age_root = sps_age_dict[i][j]
		
		# init dict of outgroups for species i
		sps_out_dict[i] = dict()
		sps_nod_dict[i] = dict()
		for m,j in enumerate(sps_list):
			if n != m:
				e = phs.get_common_ancestor(i,j)   # define species subtree between species i and j
				l = e.get_leaf_names()             # which species are present in subtree e?
				d = dict()                         # modified dict of inter-species ages: 
				for k in sps_list:                 #
					if np.isin(k, l):              # if species in subtree, use real age
						d[k] = sps_age_dict[i][k]  #
					else:                          # if species not in subtree, age = 0 
						d[k] = 0
				o = e.get_farthest_oldest_leaf(d).name # find farthest oldest leaf by name (outgroup to sps i)
				sps_out_dict[i][j] = i+","+o           # store age as pair that defines a last common ancestor
				sps_nod_dict[i][j] = phs.get_common_ancestor(i,j).name
			if n == m:
				sps_out_dict[i][j] = i+","+i      
				sps_nod_dict[i][j] = i

	return sps_age_dict, sps_out_dict, sps_nod_dict, sps_list, age_root


# loop through list of orthogroups and calculate ages
def calculate_origin(ort, sps_out_dict, sps_nod_dict, clus_col=clus_col, gene_col=gene_col, split_ch=split_ch):

	# list of orthoclusters
	clus_lis = np.unique(ort[clus_col])

	# empty lists for outputs (ages and named ages)
	ages_lis = np.zeros(len(clus_lis))
	outg_lis = np.zeros(len(clus_lis), dtype=object)
	node_lis = np.zeros(len(clus_lis), dtype=object)

	# absolute origin
	dat = pd.DataFrame()
	for n,c in enumerate(clus_lis):
		nod_clu = ort[ort[clus_col] == c][gene_col].values              # genes present in cluster
		sps_clu = np.unique([ m.split(split_ch)[0] for m in nod_clu ])  # species present in cluster

		a = 0 # init age: 0 (most recent)
		for si in sps_clu:
			for sj in sps_clu:
				t  = sps_age_dict[si][sj] # divergence time of current species pair
				if t >= a:                 # if divergence time of current pair is older than age, reassign age (a) and named ages (o)
					a = int(t)
					oi = si
					oj = sj

		ages_lis[n] = a
		outg_lis[n] = sps_out_dict[oi][oj]
		node_lis[n] = sps_nod_dict[oi][oj]

	# prepare output
	dat = pd.DataFrame({
		"orthogroup"   : clus_lis,
		"relative_age" : ages_lis,
		"outgroup_sps" : outg_lis,
		"node_age"     : node_lis
	}, 
	columns=["orthogroup", "node_age", "outgroup_sps", "relative_age" ])

	return dat


#### MAIN WORK ####

# load input tree
print("# Load tree from %s" % phs_fn)
phs    = ete3.PhyloTree("%s" % (phs_fn), format=1)
# assign species names to tree
phs.set_species_naming_function(lambda node: node.name)

# obtain species-to-species dictionary of relative ages
print("# Species-to-species relative ages from %s" % phs_fn)
sps_age_dict, sps_out_dict, sps_nod_dict, sps_list, age_root = species_age_dict(phs=phs)

# load orthoclusters
print("# Load orthoclusters from %s" % ort_fn)
ort = pd.read_csv(ort_fn, sep="\t")
ort = ort[[gene_col,clus_col]]

print("# Calculate origins...")
dat = calculate_origin(ort=ort, sps_out_dict=sps_out_dict, sps_nod_dict=sps_nod_dict)

# save output: per orthogroup
print("# Save output")
dat.to_csv("%s.node_ages_per_orthogroup.csv" % out_fn, sep="\t", index=None, mode="w")

# save output: per gene
dag = ort.merge(dat, left_on=clus_col, right_on="orthogroup")
dag = dag.drop(labels=clus_col, axis=1)
dag.to_csv("%s.node_ages_per_gene.csv" % out_fn, sep="\t", index=None, mode="w")

print("# Fet!")


