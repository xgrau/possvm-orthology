# TALE homeobox classification

1. Same species as Orthobench

2. Get TALE (ANTP) sequences from HomeoDB.

3. Blast (diamond).

4. Tree.

5. POSSVM.

```bash
python ../../possvm.py -i results_trees/tale.genes.iqtree.treefile -p tale.possom -outgroup outgroups.txt
```

6. Evaluate using classification from blast to HomeoDB.
