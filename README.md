# Automatic orthology

Scripts for automatic orthology assignment, proves.

## Pipeline

Steps:

0. Create phylogeny somehow.

1. **`findog_s01_ETEnet_v03_27set19.py`** (Python 3.5.5): takes as input a newick phylogeny and identifies speciation events using `ete`. Then, it creates a network-like table of phylogeny nodes (leaves) that appeared via the same speciation event (using various species overlap score thresholds). Each table entry is a pair of `in_seqs` and `out_seqs`, like this:

```bash
seq1	seq2	S	100.0	0.01
seq1	seq3	S	100.0	0.01
...
```

```bash
findog_s01_ETEnet_v03_27set19.py <input newick> <output prefix>
```

2. **`findog_s02_tabulate_v01.R`** (R 3.6.1): Assign each component of this network of orthologs to one or more orthogroups. Requires `igraph`.

Required inputs:

* `set_raxml.newick`: phylogeny, newick, includes supports.
* `set_Hsap_names.dict`: dictionary linking specific sequences from one specie (or more?) to gene names (orthogroups that'll be tabulated)
* `TODO`: list of species to create report (so as to include also species that are not present in the phylogeny).

Final outputs (maybe):

* table of ortholog presence/absence per species, or table of counts (including inparalogs?)
* table of ortholog support per species (branch support of the paralog with the highest support per species?)
* lists of genes belonging to each cluster

## ETE3

Two strategies available [here](http://etetoolkit.org/docs/latest/tutorial/tutorial_phylogeny.html#detecting-evolutionary-events).

### Species overlap

See [Huerta-Cepas 2007](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-6-r109). 

* Does not require a species tree.

* Only one tweakable parameter, **species overlap score** (`sos_thr`). By default, `sos_thr=0` which means that a single species in common between two node branches will rise a duplication event. This has been shown to perform the best with real data ([Huerta-Cepas 2007](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-6-r109)).

Script `ete-proves/test_v02_speciesoverlap.py` implements this method with a tree of ADAR enzymes (multigene, paralogy in animals, few outgroups). Observations:

* All identified events make sense (paralogs, orthologs) for ADAR.
* False negatives?

```bash
D Anocul_ACUA003793-RA_2-103,Anomin_AMIN009728-RA_330-722 <====> Anocul_ACUA021577-RA_55-316
```

### Strict tree reconciliation

Too strict. It'll force the species tree onto every possible subtree and add putative losses; which looks great in theory but requires a degree of compliance between gene trees and species trees that is unreasonable even for single-copy gene families. Won't for large multigene families.

## UPHO

### Setup

1. Download UPHO:

```bash
git clone https://github.com/ballesterus/UPhO.git
```

2. Run UPHO:

```bash
python2 /home/xavi/Programes/UPhO/UPhO.py -d _ -in tree/adar.newick -iP -ouT
```

Parameters:

* `-iP` to include inparalogs in the orthplogy groups
* `-m` to specify the minimum number of OTUs in an orthogroup
* `-ouT` Write orthologous branches to newick file
* `-S` minimum support value for orthology evaluation
* `-d` delimiter, default is `|`, but if you use `_` it seems to work fine (even if there are more than one such characters per name)

Problem: highly atomised for some reason! There's no way to tweak granularity.

## GeneRax

Check [this](https://www.biorxiv.org/content/10.1101/779066v1).

## Create orthogroups, alignments and trees

First, run orthofinder:

```bash
orthofinder -t 4 -a 4 -M msa -S diamond -A mafft -T iqtree -f input/ -I 1.5 -s tree.newick -os    ### -os ensures that sequence files are created
```

Create long-format orthofinder output:

```bash
awk '{ for (i=2; i <= NF; i++) { print $i"\t"$1  }}' Orthogroups_nomllarg.txt | sed "s/://" > Orthogroups_nomllarg_llarg.csv
```


Create `Alignments` and `Trees` folders:

```bash
cd input/Orthologes_XX/Sequences
mkdir ../Alignments
mkdir ../Trees
```

Run alignments and trees:

* Checks if orthogroups contain >1 seq (otherwise, alignment fails)
* Runs mafft, trimal and fasttree (iqtree is too slow?)

```bash
for i in *.fa ; do if [ $(grep -c ">" $i) -gt 1 ] ; then if [ ! -f ../Alignments/${i%%.fa}.l.fa ]; then echo ${i%%.fa} ali ; mafft --localpair --reorder --maxiterate 1000 --thread 6 $i > ../Alignments/${i%%.fa}.l.fa 2> /dev/null ; fi ; if [ ! -f ../Alignments/${i%%.fa}.lt.fa ] ; then trimal -in ../Alignments/${i%%.fa}.l.fa -out ../Alignments/${i%%.fa}.lt.fa -automated1 ; fi ; if [ ! -f ../Fasttrees/${i%%.fa}.tree ] ; then echo ${i%%.fa} phy ;  fasttree -lg -quiet -cat 4 ../Alignments/${i%%.fa}.lt.fa > ../Fasttrees/${i%%.fa}.tree ; fi ; fi ; done
```

Then, find a way to combine this output with `ETE` + species overlap.

* Currently, script in `orthofinder_Ano14sp`.
* Runs `ETE` + species overlap and finds clusters with `MCL` (inflation 1.5).

```bash
python findog_s01_ETEmcl_v07_10oct19.py /home/xavi/dades/Anotacions/orthofinder_Ano14sps_noclu_9oct19/output/Fasttrees out10oct19
```

Obtain species tree with distances, using single-copy orthologs:

```bash
while read s; do echo ">$s" ; for i in $(cat SingleCopyOrthogroups.txt) ; do sed "s/\(>[^_]\)*_.*/\1/" Alignments/${i}.lt.fa | bioawk -c fastx '{ print $1,$2 }' | grep "^$s" | cut -f2 ; done   ; done < ../tree.list  > ../tree_dists.ali.fasta ; iqtree -g ../tree.newick.unroot -m LG+G4 -s ../tree_dists.ali.fasta -pre ../tree_dists.iqt -redo
```

Now, parse with `phytools`: <http://www.phytools.org/eqg2015/asr.html>
