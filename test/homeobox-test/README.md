# TALE homeobox classification

1. Get TALE, PRD and ANTP sequences from [HomeoDB](http://homeodb.zoo.ox.ac.uk/) (as of 7th Feb 2021).

2. Blast (diamond).

```bash
# to orthobench dataset (17 metazoan species)
bash s01_get_trees-diamond.sh seed_tale.fasta tale ../orthobench-test/proteomes/
bash s01_get_trees-diamond.sh seed_ANTP.fasta ANTP ../orthobench-test/proteomes/
bash s01_get_trees-diamond.sh seed_PRD.fasta PRD ../orthobench-test/proteomes/
# to Stylophora dataset (35 metazoan species)
bash s01_get_trees-diamond.sh seed_ANTP.fasta ANTPm proteomes/
bash s01_get_trees-diamond.sh seed_PRD.fasta PRDm  proteomes/
bash s01_get_trees-diamond.sh seed_tale.fasta TALEm  proteomes/
```

3. POSSVM.

```bash
python ../../possvm.py -i results_trees/tale.genes.iqtree.treefile -p tale.possom -outgroup outgroups.txt
```

4. Evaluate using classification from blast to HomeoDB.

```bash
Rscript s02_evaluate_homeodb.R
```

5. Plot and annotate trees?
