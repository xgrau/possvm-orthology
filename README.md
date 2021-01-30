# Possvm

**Possvm** (***P**hylogenetic **O**rtholog **S**orting with **S**pecies O**v**erlap and **M**CL*) is a python utility that analyses pre-computed gene trees to identify orthologous sequences. It takes advantage of the **[ETE toolkit](http://etetoolkit.org/)** to parse the phylogeny and identify orthologous gene pairs, and **[MCL clustering](https://micans.org/mcl/)** for orthogroup identification.

Its basic functionality only requires a gene tree in newick format, with sequence name containing a prefix that indicates their species of origin, e.g. `human_gene1`. It does *not* require a species tree to infer orthologs, because it relies on the **[species overlap algorithm](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-6-r109)** implemented in ETE (see [here](http://etetoolkit.org/docs/latest/tutorial/tutorial_phylogeny.html#species-overlap-so-algorithm)).

## How it works

Dibuix.

## How to cite

Please cite the following papers:

* *POSSVM* paper: **[here]**.
* *ETE* toolkit: **[Huerta-Cepas *et al.* Molecular Biology and Evolution 2016](http://etetoolkit.org/)**.
* Species overlap algorithm: **[Huerta-Cepas *et al.* Genome Biolgy 2007](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-6-r109)**.
* *MCL* clustering: **[Enright *et al.* Nucleic Acids Research 2002](https://micans.org/mcl/)**.

## Manual

### Main script: `possvm.py`

**`possvm.py`** is the main script, used to parse phylogenetic trees and obtain sub-groups of genes that constitute ortholog clusters.

Usage:

```bash
usage: possvm.py [-h] -p PHY -o OUT [-i ID] [-r REF]
                 [-refsps REFSPS] [-s SOS]
                 [-split SPLIT] [-skiproot]
                 [-skipprint]
                 [-min_transfer_support MIN_TRANSFER_SUPPORT]
                 [-extratio EXTRATIO]

optional arguments:
  -h, --help            show this help message and exit
  -p PHY, --phy PHY     Path to a phylogenetic tree in newick format. Each
                        sequence in the tree must have a prefix indicating the
                        species, separated from gene name with a split
                        character. Default split character is "_", see --split
                        for options.
  -o OUT, --out OUT     Path to output folder. Defaults to present working
                        directory.
  -i ID, --id ID        OPTIONAL: String. Gene family name, used when naming
                        ortholog clusters. Defaults to "genefam".
  -r REF, --ref REF     OPTIONAL: Path to a table indicating reference gene
                        names that can be used for orthogroup labeling.
                        Format: geneid <tab> name.
  -refsps REFSPS, --refsps REFSPS
                        OPTIONAL: Comma-separated list of reference species
                        that will be used for orthogroup labeling.
  -s SOS, --sos SOS     OPTIONAL: Species overlap threshold used for orthology
                        inference in ETE. Default is 0.
  -split SPLIT, --split SPLIT
                        OPTIONAL: String to use as species prefix delimiter in
                        gene ids, e.g. "_" for sequences formatted as
                        speciesA_geneX. Defaults to "_".
  -skiproot, --skiproot
                        OPTIONAL: Turns off tree rooting using midpoint root,
                        in case your trees are already rooted.
  -skipprint, --skipprint
                        OPTIONAL: Turns off printing of annotated tree in PDF
                        (annotated newick is still produced).
  -min_transfer_support MIN_TRANSFER_SUPPORT, --min_transfer_support MIN_TRANSFER_SUPPORT
                        OPTIONAL: Min node support to allow transfer of labels
                        from labelled to non-labelled groups in the same
                        clade. If not set, this step is skipped.
  -extratio EXTRATIO, --extratio EXTRATIO
                        NOT IN USE!! OPTIONAL: In order to perform extended
                        label propagation, you can assign XX. Ratio Defaults
                        to 1.5, ie closest group is 50pp loser to unlabelled
                        group than the second closest group.

```

##### Input

Phylogenies must be in **newick format** and can contain node supports.

If you are using the **single tree mode**, use `-phy` to point to the gene tree.

If you are using the **tree collection mode**, use `-phy` to point to the tree folder. Each phylogeny in the tree folder must be named as follows: `orthogroup_name.suffix`. The suffix is indicated with the `-suf` flag. Examples:

* `OG00001.newick`: `-suf newick`
* `OG00001.iqtree.treefile`: `-suf iqtree.treefile`

#### Output

* describe

### Gene ages with `possvm_geneage.py`

**`possvm_geneage.py`**: a simple way to obtain the **ages of genes and clusters of orthologs**. It takes as input the output table from `possvm.py` (or any similarly formatted table) and a species tree, and outputs a new table with the age of each orthogroup and each gene.

All gene ages are relative to a pre-defined species of reference, and are defined as follows:

* relative numeric
* relative outgroup pairs
* named ancestors, if available

It can also be used to analyse outputs from Orthofinder/OrthoMCL. 

This script **requires a species tree** in newick format.

Usage:

```bash
python possvm_nodeage.py -h
usage: possvm_nodeage.py [-h] -tree TREE -ort ORT -out OUT [-ref REF]
                         [-dict DICT] [-gcol GCOL] [-ccol CCOL] [-split SPLIT]

optional arguments:
  -h, --help            show this help message and exit
  -tree TREE, --tree TREE
                        newick file with species tree
  -ort ORT, --ort ORT   table of orthologs (long format, with headers; min
                        info: OG <tab> gene1)
  -out OUT, --out OUT   output prefix
  -ref REF, --ref REF   OPTIONAL: reference species (ages are relative to this
                        species). If "ALLSPS" (default), relative ages to all
                        sps are calculated (can be slow)
  -dict DICT, --dict DICT
                        OPTIONAL: dictionary file with named nodes, one per
                        line (can be incomplete). Format: spsA,spsB <tab>
                        nodename. If "NO" (default), do nothing
  -gcol GCOL, --gcol GCOL
                        OPTIONAL: name of gene column in --ort/-t table.
                        Default is "node"
  -ccol CCOL, --ccol CCOL
                        OPTIONAL: name of orthogroup column in --ort/-t table.
                        Default is "og_cluster"
  -split SPLIT, --split SPLIT
                        OPTIONAL: character to split species and sequence
                        names. Default is "_", e.g. Human_genename. WARNING:
                        use quotation marks, e.g. -split "_" or -split "|"

```

### Examples

#### Analyse a collection of gene trees

Run **`possvm.py`** on a collection of trees:

```bash
python possvm.py -phy test_anopheles/trees/ -suf newick -out test_anopheles/output_etemcl -ort test_anopheles/Orthogroups_longformat.csv -ani main
```

Check output in `test_anopheles/output_etemcl.orthology.csv`:

```bash
$ head test_anopheles/output_etemcl.orthology.csv
node	og	cluster	og_cluster
Anoalb_AALB001334-RA	OG0000000		OG0000000_nan
Anoalb_AALB002504-RA	OG0000000	10	OG0000000_10
Anoalb_AALB002511-RA	OG0000000	2	OG0000000_2
Anoalb_AALB002512-RA	OG0000000	3	OG0000000_3
Anoalb_AALB003039-RA	OG0000000	4	OG0000000_4
Anoalb_AALB008232-RA	OG0000000	9	OG0000000_9
Anoalb_AALB008239-RA	OG0000000	12	OG0000000_12
Anoalb_AALB008240-RA	OG0000000	20	OG0000000_20
Anoalb_AALB015485-RA	OG0000000	0	OG0000000_0
```

#### Analyse a a single gene tree

Run **`possvm.py`** on a collection of trees:

```bash
python possvm.py -phy test_anopheles/trees/ -suf newick -out test_anopheles/output_etemcl -ort test_anopheles/Orthogroups_longformat.csv -ani main
```

Check output in `test_anopheles/output_etemcl.orthology.csv`:

```bash
$ head test_anopheles/output_etemcl.orthology.csv
node	og	cluster	og_cluster
Anoalb_AALB001334-RA	OG0000000		OG0000000_nan
Anoalb_AALB002504-RA	OG0000000	10	OG0000000_10
Anoalb_AALB002511-RA	OG0000000	2	OG0000000_2
Anoalb_AALB002512-RA	OG0000000	3	OG0000000_3
Anoalb_AALB003039-RA	OG0000000	4	OG0000000_4
Anoalb_AALB008232-RA	OG0000000	9	OG0000000_9
Anoalb_AALB008239-RA	OG0000000	12	OG0000000_12
Anoalb_AALB008240-RA	OG0000000	20	OG0000000_20
Anoalb_AALB015485-RA	OG0000000	0	OG0000000_0
```

### Requirements

Also:
