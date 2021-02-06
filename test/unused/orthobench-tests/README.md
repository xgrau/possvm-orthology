# Orthobench tests

1. Get Orthobench reference groups from the [Open_Orthobench Github](https://github.com/davidemms/Open_Orthobench/tree/master/BENCHMARKS) (relevant files in `BENCHMARKS/`: `RefOGs`, `Input` and `benchmark.py`; total: 70 RefOGs).

2. NOTE: Fix duplicated assignment of `FBpp0309618` to the `RefOG021` and `RefOG068` clusters (remove from latter).

3. Get Orthobench reference gene trees from [Open_Orthobench Github](https://github.com/davidemms/Open_Orthobench/tree/master/Supporting_Data/Data_for_RefOGs/trees); remove provisional versions (total: 70 trees).

4. Run *POSSVM* as follows:

```bash
for i in tree-collection/raw/*.tre ; do python ../../possom.py -i $i -o tree-collection-parsed -skiproot -ogprefix "$(basename $i | sed "s/.tre//")." ; done
```

5. Prepare table for benchmark:

```bash
# create table with all orthogroups from trees
for i in tree-collection/raw-parsed/*_groups.csv ; do awk 'NR >1' $i | sed -r "s/^[A-Za-z]+_[A-Za-z]+_//" | awk '{ if(a[$2]) a[$2]=a[$2]" "$1 ; else a[$2]=$1 ; }END { for (i in a) { print i":\t" a[i]; } }' ; done | shuf  > possom_RefOGs.csv
# (trees contain redundant genes, so I'll sample random first hits to remove redundancy and retain orthogroups that contain genes in the reference dataset)
# for i in RefOGs/*.txt ; do fgrep -w -f $i -m 1 possom_RefOGs.csv ; done | sort -u > possom_RefOGs_unique.csv
# alternative: pick group that contains the highest number of REFOG genes:
for i in RefOGs/*.txt ; do fgrep -n -w -f $i possom_RefOGs.csv | sed "s/:/\t/"  | sort -k1,1nr | head -n1 | cut -f2,3 ; done | sort -u > possom_RefOGs_unique.csv
# run benchmark
python benchmark.py possom_RefOGs_unique.csv

```

6. Benchmark output:

```
1780 genes found in predicted orthogroups, this is 91.6% of all genes

Calculating benchmarks:
76.7% F-score
79.1% Precision
74.4% Recall

43 Orthogroups exactly correct
```

DOES ANY OF THIS MAKE SENSE?

## Alternative

Build trees from input dataset, with OrthoFinder:

```bash
python ~/Programes/OrthoFinder/orthofinder.py -t 8 -f Input/ -M msa -T iqtree -A mafft -o output_mafft_iqtree
python ~/Programes/OrthoFinder/orthofinder.py -t 8 -f Input/ -M dendroblast -o output_dendroblast
```

And then parse these trees with *POSSVM*, and calculate benchmarks as above.
