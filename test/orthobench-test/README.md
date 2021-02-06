# Orthobench testing

Test the accuracy of *Possvm* using manually curated orthologs from the [Orthobench 2.0 repository](https://github.com/davidemms/Open_Orthobench) (see [Emms et al. GBE 2020](https://academic.oup.com/gbe/article/12/12/2258/5918455)).

# Steps

1. Get Orthobench reference groups (`refOG`) from the [Open_Orthobench Github](https://github.com/davidemms/Open_Orthobench/tree/master/BENCHMARKS); total: 70 RefOGs), as well as raw and tight trees (ignore intermediate trees).
:

```bash
# get Orthobench repository
git clone git@github.com:davidemms/Open_Orthobench.git

# refOGs assignments:
ls <path>/Open_Orthobench/BENCHMARKS/RefOGs/

# tight and raw trees:
ls <path>/Open_Orthobench/Supporting_Data/Data_for_RefOGs/trees
ls <path>/Open_Orthobench/Supporting_Data/Data_for_RefOGs/trees_tight

# # HMM profiles in this folder
# ls <path>/Open_Orthobench/Supporting_Data/Additional_Files/hmm_profiles/
# # copy the HMM files into the `hmm_profiles_strict` and `hmm_profiles_weak` folders in the present directory
```

2. NOTE: Fix duplicated assignment of `FBpp0309618` to the `RefOG021` and `RefOG068` clusters (remove from latter), format into `refOGs.csv`.

3. Get primary transcripts from 17 metazoan species from the [Orthobench repository](https://github.com/davidemms/Open_Orthobench/tree/master/Supporting_Data/Additional_Files/proteomes/primary_transcripts), and HMM profiles for each refOG:

```bash
# primary transcripts for each species, in this folder:
ls <path>/Open_Orthobench/Supporting_Data/Additional_Files/proteomes/primary_transcripts
# copy the fasta files into the `proteomes` folder in in the present directory and concatenate them...
```

4. Run homology searches, MSAs and phylogenies:

```bash
# run from the present directory
bash s01_get_trees-diamond.sh
```

5. Run *Possvm* as follows:

```bash
# find OGs in each tree with possom:
for i in orthobench_trees/tight/*.tre ; do python ../../possvm.py -i $i -p  $(basename ${i%%.*}).possom -ogprefix "$(basename ${i%%.*})." -outgroup outgroups.txt ; done
for i in orthobench_trees/raw/*.tre ; do python ../../possvm.py -i $i -p  $(basename ${i%%.*}).possom -ogprefix "$(basename ${i%%.*})." -outgroup outgroups.txt ; done
for i in results_orthology/*.treefile ; do python ../../possvm.py -i $i -p  $(basename ${i%%.*}).possom -ogprefix "$(basename ${i%%.*})." -outgroup outgroups.txt ; done
```

6. Calculate accuracy relative to `refOGs.csv`:

```bash
Rscript s02_evaluate_orthobench.R
```
