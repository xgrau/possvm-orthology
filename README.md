# Possvm

**Possvm** (_**P**hylogenetic **O**rtholog **S**orting with **S**pecies O**v**erlap and **M**CL_) is a python utility that analyses pre-computed gene trees to identify orthologous sequences. It takes advantage of the **[*ETE* toolkit](http://etetoolkit.org/)** to parse the phylogeny and identify orthologous gene pairs, and **[*MCL* clustering](https://micans.org/mcl/)** for orthogroup identification.

Its basic functionality only requires a gene tree in newick format, with sequence name containing a prefix that indicates their species of origin, e.g. `human_gene1`. It does *not* require a species tree to infer orthologs, because it relies on the **[species overlap algorithm](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-6-r109)** implemented in ETE (see [here](http://etetoolkit.org/docs/latest/tutorial/tutorial_phylogeny.html#species-overlap-so-algorithm)).

## How it works

Dibuix.

## How to cite

If you use *Possvm*, please cite the following papers:

* *Possvm* paper: **[here]**.
* *ETE* toolkit: **[Huerta-Cepas *et al.* Molecular Biology and Evolution 2016](http://etetoolkit.org/)**.
* Species overlap algorithm: **[Huerta-Cepas *et al.* Genome Biolgy 2007](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-6-r109)**.
* *MCL* clustering: **[Enright *et al.* Nucleic Acids Research 2002](https://micans.org/mcl/)**.

## Manual

### Main script: `possvm.py`

**`possvm.py`** is the main script, used to parse phylogenetic trees and obtain sub-groups of genes that constitute ortholog clusters.

Usage:

```man
usage: possvm.py [-h] -i IN [-o OUT] [-p PHY] [-r REF] [-refsps REFSPS]
                 [-s SOS] [-outgroup OUTGROUP] [-split SPLIT]
                 [-itermidroot ITERMIDROOT] [-skiproot] [-skipprint]
                 [-min_transfer_support MIN_TRANSFER_SUPPORT]
                 [-clean_gene_names] [-cut_gene_names CUT_GENE_NAMES]
                 [-ogprefix OGPREFIX]

optional arguments:
  -h, --help            show this help message and exit
  -i IN, --in IN        Path to a phylogenetic tree in newick format. Each
                        sequence in the tree must have a prefix indicating the
                        species, separated from gene name with a split
                        character. Default split character is "_", see --split
                        for options.
  -o OUT, --out OUT     OPTIONAL: Path to output folder. Defaults to same
                        directory as input file.
  -p PHY, --phy PHY     OPTIONAL: Prefix for output files. Defaults to
                        `basename` of input phylogeny. Default behaviour will
                        never overwrite original files, because it adds
                        suffixes.
  -r REF, --ref REF     OPTIONAL: Path to a table indicating reference gene
                        names that can be used for orthogroup labeling.
                        Format: geneid <tab> name.
  -refsps REFSPS, --refsps REFSPS
                        OPTIONAL: Comma-separated list of reference species
                        that will be used for orthogroup labeling. If absent,
                        all sequences present in the -r table will be
                        considered.
  -s SOS, --sos SOS     OPTIONAL: Species overlap threshold used for orthology
                        inference in ETE. Default is 0.
  -outgroup OUTGROUP, --outgroup OUTGROUP
                        OPTIONAL: String. Define a set of species that are
                        treated as outgroups in the phylogeny, and excluded
                        from orthology clustering. Can be a comma-separated
                        list of species, or a file with one species per line.
                        This option DOES NOT affect tree rooting, just
                        orthology clustering. Disabled by default.
  -split SPLIT, --split SPLIT
                        OPTIONAL: String to use as species prefix delimiter in
                        gene ids, e.g. "_" for sequences formatted as
                        speciesA_geneX. Defaults to "_".
  -itermidroot ITERMIDROOT, --itermidroot ITERMIDROOT
                        OPTIONAL: Turns on iterative midpoint rooting with N
                        iterations, which is used instead of the default
                        midpoint rooting.
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
  -clean_gene_names, --clean_gene_names
                        OPTIONAL: Will attempt to "clean" gene names from the
                        reference table (see -r) used to create cluster names,
                        to avoid very long strings in groups with many
                        paralogs. Currently, it collapses number suffixes in
                        gene names, and converts strings such as Hox2/Hox4 to
                        Hox2-4. More complex substitutions are not supported.
  -cut_gene_names CUT_GENE_NAMES, --cut_gene_names CUT_GENE_NAMES
                        OPTIONAL: Integer. If set, will shorten cluster name
                        strings to the given length in the PDF file, to avoid
                        long strings in groups with many paralogs. Default is
                        no shortening.
  -ogprefix OGPREFIX, --ogprefix OGPREFIX
                        OPTIONAL: String. Prefix for ortholog clusters.
                        Defaults to "OG".
```

Input:

* Phylogenies must be in **newick format** and can contain node supports.

* If you are using the **single tree mode**, use `-phy` to point to the gene tree.

* If you are using the **tree collection mode**, use `-phy` to point to the tree folder. Each phylogeny in the tree folder must be named as follows: `orthogroup_name.suffix`. The suffix is indicated with the `-suf` flag.

Output:

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


### Install *Possvm* and its dependencies

***Possvm*** depends on the [*ETE* toolkit](http://etetoolkit.org/) Python library, which currently works best with Python 3.6 and can be installed via *conda*. We thus recommend that you use conda to install *ETE* and all other dependencies.

Once you have a working installation of *conda* (see [here for instructions](http://etetoolkit.org/download/)), you can run the following commands:

```bash
# create environment for possvm
conda create -n possvm python=3.6
# install ETE toolkit
conda install -c etetoolkit ete3
# install other dependencies
conda install -c bioconda pandas networkx markov_clustering matplotlib
# activate the environment
conda activate possvm
```

Alternatively, you can use the `environment.yaml` file bundled in this repository to reproduce the environment:

```bash
# create env and install packages
conda env create -n possvm --file environment.yaml
# activate the environment
conda activate possvm
```

Both options should download and install all basic dependencies, including the following packages:

```bash
ete3                      3.1.2              pyh39e3cac_0    etetoolkit
markov_clustering         0.0.6                      py_0    bioconda
matplotlib                3.3.2                h06a4308_0  
matplotlib-base           3.3.2            py36h817c723_0  
networkx                  2.5                        py_0  
numpy                     1.19.2           py36h54aff64_0  
numpy-base                1.19.2           py36hfa32c7d_0  
pandas                    1.1.3            py36he6710b0_0  
python                    3.6.12               hcff3b4d_2  
```

Once these dependencies are up and running, you can run *Possvm* like any Python script:

```bash
python possvm.py -h

# Or maybe add it as an alias?
echo "alias possvm=\"python $(pwd)/possvm.py\"" >> ~/.bashrc
source ~/.bashrc
possvm -h
```

### Examples

Some examples

