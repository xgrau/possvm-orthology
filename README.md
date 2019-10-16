# Possom

**Possom** (***P**hylogenetic **O**rtholog **C**lustering with **S**pecies **O**verlap and **M**CL*) is a python utility that analyses gene phylogenies and defines clusters of orthologs within each tree, taking advantage of the **[ETE toolkit](http://etetoolkit.org/)** for phylogeny analysis and **[MCL clustering](https://micans.org/mcl/)**.

It can be used to analyse a single gene tree or a collection of gene trees (for example, gene trees obtained from orthogroups defined using [Orthofinder](https://github.com/davidemms/OrthoFinder)).

It requires two very basic inputs (see details below):

* One or more gene trees
* If more than one gene tree is used, the user should specify a list of 

Its basic functionality does *not* require a species tree. This is because it relies on the **[species overlap algorithm](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-6-r109)** implemented in ETE (see [here](http://etetoolkit.org/docs/latest/tutorial/tutorial_phylogeny.html#species-overlap-so-algorithm)).

## How it works

Dibuix.

## How to cite

Please cite the following papers:

* Mine.
* ETE toolkit: **[Huerta-Cepas *et al.* Molecular Biology and Evolution 2016](http://etetoolkit.org/)**.
* Species overlap algorithm: **[Huerta-Cepas *et al.* Genome Biolgy 2007](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-6-r109)**.
* MCL clustering: **[Enright *et al.* Nucleic Acids Research 2002](https://micans.org/mcl/)**.

## Manual

### Main script: `possom.py`

**`possom.py`** is the main script, used to parse phylogenetic trees and obtain sub-groups of genes that constitute ortholog clusters.

Usage:

```bash
$ python possom.py -h
usage: possom.py [-h] -phy PHY -suf SUF -out OUT -ort ORT -ani ANI [-inf INF]
                 [-nopt NOPT] [-print PRINT] [-split SPLIT]

optional arguments:
  -h, --help            show this help message and exit
  -phy PHY, --phy PHY   String. Folder with phylogenies
  -suf SUF, --suf SUF   String. Suffix of phylogenies in folder
  -out OUT, --out OUT   String. Prefix for output
  -ort ORT, --ort ORT   String. Path to orthology file. Must be a two-column
                        table (with tabs), with one gene per line: OG <tab>
                        gene1
  -ani ANI, --ani ANI   String. Which analysis to perform: "main" to analyse
                        all genes, "opti" to find optimal inflation value
  -inf INF, --inf INF   OPTIONAL: Floating. if analysis is "main", which
                        inflation value to use? Default is 1.1
  -nopt NOPT, --nopt NOPT
                        OPTIONAL: Integer. if analysis is "opti", how many
                        phylogenies should we examine for optimisation?
                        Default is 500
  -print PRINT, --print PRINT
                        OPTIONAL: Boolean (False/True). Print new tree with
                        defined clusters?
  -split SPLIT, --split SPLIT
                        OPTIONAL: character to split species and sequence
                        names. Default is "_", e.g. Human_genename. WARNING:
                        use quotation marks, e.g. -split "_" or -split "|"
```

#### Inputs

* Phylogenies in the tree folder (`-phy`) must be in newick format and can contain bootstrap supports, that will be used for MCL clustering
* Each phylogeny in the tree folder should be named as follows: `orthogroup_name.suffix`. The suffix is indicated with the `-suf` flag. For example: `OG00001.newick` (`-suf newick`) or `OG00001.iqtree.treefile` (`-suf iqtree.treefile`).
* The **table of orthologs** (`-ort`) must be formatted as follows:
```
gene1	OG1
gene2	OG1
gene3	OG2
gene4	OG2
...
```

If you obtained orthologs using Orthofinder, you can use this awk one-liner to obtain a table formatted as above from the `Orthofinder.txt` file:

```bash
awk '{ for (i=2; i <= NF; i++) { print $i"\t"$1 }}' Orthogroups.txt | sed "s/://" > Orthogroups_longformat.csv
```
* 

Examples of compliant input fles are provided in the `test_anopheles` folder.

### Gene ages: `possom_geneage.py`

**`possom_geneage.py`**: a simple way to obtain the ages of each cluster of orthologs identified in the phylogeny. It can also be used to analyse outputs from Orthofinder/OrthoMCL. 

This script **requires a species tree** in newick format.

Usage:

```bash
python possom_nodeage.py -h
/home/xavi/Programes/miniconda3/envs/pyco/lib/python3.5/importlib/_bootstrap.py:222: RuntimeWarning: numpy.dtype size changed, may indicate binary incompatibility. Expected 96, got 88
  return f(*args, **kwds)
usage: possom_nodeage.py [-h] -tree TREE -ort ORT -out OUT [-ref REF]
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

Run **`possom.py`** on a collection of trees:

```bash
python possom.py -phy test_anopheles/trees/ -suf newick -out test_anopheles/output_etemcl -ort test_anopheles/Orthogroups_longformat.csv -ani main
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

Run **`possom.py`** on a collection of trees:

```bash
python possom.py -phy test_anopheles/trees/ -suf newick -out test_anopheles/output_etemcl -ort test_anopheles/Orthogroups_longformat.csv -ani main
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

## Requirements


Also