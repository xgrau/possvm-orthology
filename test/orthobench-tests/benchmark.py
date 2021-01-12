#!/usr/bin/env python3
"""
This script calculates the benchmarks for an input set of orthogroups

Instructions:
1. Predict the complete set of orthogroups for the genes in the "Input/" directory
2. Write Orthogroups to a file, one orthogroup per line (header lines and commented
   out lines starting with '#' are allowed). 
3. Download this script and the accompanying RefOGs directory
3. Call the script with the the orthogroup filename as the only argument

By default the script with use regular expressions to extract the genes from the 
additional text on each line. 

You can also specify the option '-b' to use the more basic file reader which requires 
the following format:
- First line is a header and is ignored. 
- One orthogroup per line
- Genes can be separated by commas, spaces or tabs. 
- Lines starting with '#' are comments and are ignored
- Optionally, each line can start with the name of the orthogroup followed by a colon.
   E.g. "OG0000001: Gene1, Gene2"
"""
import os
import re
import sys
import glob
import argparse


delims = " |,|\t"
expected_gene_names_base = {"ENSP000": "Homo sapiens", 
    "ENSRNO":"Rattus norvegicus", 
    "ENSCAF":"Canis familiaris", 
    "ENSMUS":"Mus musculus", 
    "ENSMOD":"Monodelphis domestica", 
    "ENSGAL":"Gallus gallus", 
    "ENSCIN":"Ciona intestinalis",
    "ENSTNI":"Tetraodon nigroviridis",
    "ENSPTR":"Pan troglodytes",
    "ENSDAR":"Danio rerio",
    "FBpp":"Drosophila melanogaster",
    "WBGene":"C elegans"}
n_genes_total = 1944


def exit():
    sys.exit()


def read_hierarchical_orthogroup(infile):
    ogs = []
    iSkip = 3 
    for line in infile:
        genes = [g for s in line.rstrip().split("\t")[iSkip:] for g in s.split(", ") if g != ""]
        ogs.append(set(genes))
    return ogs

def read_orthogroups_smart(fn):
    ogs = []
    gene_pat = re.compile("WBGene00\d+\.1|ENSCAFP\d+|ENSCINP\d+|ENSDARP\d+|FBpp0\d+|ENSGALP\d+|ENSP000\d+|ENSMODP\d+|ENSMUSP\d+|ENSPTRP\d+|ENSRNOP\d+|ENSTNIP\d+")
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith("#"):
                continue
            genes = re.findall(gene_pat, l)
            # print(genes)
            # sys.exit()
            if len(genes) > 0:
                ogs.append(set(genes))
    return ogs

def read_orthogroups(fn, exp_genes, n_col_skip=1):
    """
    Read the orthogroups from a file formatted as specified above
    """
    ogs = []
    q_past_header = False
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith("#"):
                continue
            t = l.rstrip().split(None, n_col_skip)[-1]
            genes = re.split(delims, t)
            # remove white spaces
            genes = [g.strip() for g in genes]
            genes = set([g for g in genes if g != ""])
            if q_past_header or (exp_genes & genes):
                q_past_header = True
                if len(genes) > 1: ogs.append(genes)
    return ogs


def get_n_col_skip(fn):
    with open(fn, 'r') as infile:
        header = next(infile)
        if header.startswith("HOG\tOG\tGene Tree Parent Clade"):
            return 3
    return 1


def get_expected_genes():
    d_input = "Input" + os.sep
    input_files = list(glob.glob(d_input + "*fa"))
    assert(12 == len(input_files))
    all_genes = []
    for fn in input_files:
        with open(fn, 'r') as infile:
            for l in infile:
                if l.startswith(">"):
                    all_genes.append(l[1:].rstrip())
    all_genes = set(all_genes)
    assert(len(all_genes) == n_genes_total)
    return all_genes


def check_orthogroups(ogs, exp_genes):
    all_pred_genes = set([g for og in ogs for g in og])
    x = all_pred_genes.difference(exp_genes)
    if len(x) != 0:
        print("ERROR: found extra genes in input file, check its formatting is correct and there are no incorrect genes")
        print("Examples:")
        for g in list(x)[:10]:
            print(g)
    x = exp_genes.difference(all_pred_genes)
    if len(x) != 0:
        print("Examples of genes not in file:")
        for g in list(x)[:3]:
            print(g)

    n_genes = sum([len(og) for og in ogs])
    if n_genes < 0.5 * n_genes_total:
        print("ERROR: Too many missing genes in predicted orthogroups.")
        print("Orthogroups should contain at least 50% of all genes but")
        print("orthogroups file only contained %d genes" % n_genes)
        exit()
    all_genes = [g for og in ogs for g in og]
    n_genes_no_dups = len(set(all_genes))
    if n_genes_no_dups != n_genes:
        print("ERROR: Some genes appear in multiple orthogroups, benchmark are meaningless with such data.")
        print("with such data")
        print((n_genes_no_dups, n_genes))
        from collections import Counter
        c = Counter(all_genes)
        for g, n in c.most_common(1000):
            if n > 1:
                print("%d: %s" % (n,g))
        raise Exception()
    # checked genes from each of the expected species are present
    for g_pat, sp in expected_gene_names_base.items():
        if not any(g.startswith(g_pat) for g in all_genes):
            print("ERROR: No genes found from %s" % sp)
    p = 100.*n_genes/float(n_genes_total)
    print("%d genes found in predicted orthogroups, this is %0.1f%% of all genes" % (n_genes, p))



def read_refogs(d_refogs):
    refogs = []
    for i in range(1, 71):
        fn = d_refogs + ("RefOG%03d.txt" % i)
        if not os.path.exists(fn):
            print("ERROR: RefOG file not found: %s" % fn)
            exit()
        with open(fn, 'r') as infile:
            refogs.append(set([g.rstrip() for g in infile.readlines()]))
    n = sum([len(r) for r in refogs])
    n_expected = 1945
    if not n_expected == n:
        print("ERROR: There are genes missing from the RefOG files. Found %d, expected %d" % (n, n_expected))
    return refogs


def read_uncertain_refogs(d_refogs):
    refogs = []
    for i in range(1, 71):
        fn = d_refogs + ("RefOG%03d.txt" % i)
        if not os.path.exists(fn): 
            refogs.append(set())
        else:
            with open(fn, 'r') as infile:
                refogs.append(set([g.rstrip() for g in infile.readlines()]))
    return refogs


def calculate_benchmarks_pairwise(ref_ogs, uncert_genes, pred_ogs, q_even=True, q_remove_uncertain=True):
    referenceOGs = ref_ogs
    predictedOGs = pred_ogs
    totalFP = 0.
    totalFN = 0.
    totalTP = 0.
    totalGroundTruth = 0.
    n_exact = 0
    n_splits = []
    # so as not to count uncertain genes either way remove them from the 
    # expected and remove them from any predicted OG (as though they never existed!)
    for refOg, uncert in zip(referenceOGs, uncert_genes):
        thisFP = 0.
        thisFN = 0.
        thisTP = 0.
        this_split = 0
        if q_remove_uncertain:
            refOg = refOg.difference(uncert)
        nRefOG = len(refOg)
        not_present = set(refOg)
        for predOg in predictedOGs:
            overlap = len(refOg.intersection(predOg))
            if overlap > 0:
                if q_remove_uncertain:
                    predOg = predOg.difference(uncert)   # I.e. only discount genes that are uncertain w.r.t. this RefOG
                overlap = len(refOg.intersection(predOg))
            if overlap > 0:
                this_split += 1
                not_present = not_present.difference(predOg)
                thisTP += overlap * (overlap - 1)/2    # n-Ch-2
                thisFP += overlap * (len(predOg) - overlap)
                thisFN += (nRefOG - overlap) * overlap
        # Are FNs more from splintered OGs or missing genes?
        # print("%f\t%f" % (thisFN/2./(nRefOG-1), len(not_present)*(nRefOG-1)/2./(nRefOG-1))) 
        # finally, count all the FN pairs from those not in any predicted OG
        thisFN += len(not_present)*(nRefOG-1)
        # don't count 'orphan genes' as splits, it's more informative only to count 
        # clusters that this orthogroup has been split into. Recall already counts
        #  the missing genes, this should give a distinct measure.
        # this_split += len(not_present)     
        # All FN have been counted twice
        assert(thisFN % 2 == 0)
        n_splits.append(this_split)
        # print(this_split)      # Orthogroup fragments
        # print(len(not_present))  # Unclustered genes 
        thisFN /= 2 
        # sanity check
        nPairs1 = thisTP + thisFN
        nPairs2 = nRefOG * (nRefOG - 1) / 2
        if nPairs1 != nPairs2:
            print("ERROR: %d != %d" % (nPairs1,nPairs2))    
            # print(refOg)
            # print(predOg)
            # print((nRefOG, thisTP, thisFP, thisFN))
            # print("ERROR: Sanity check failed")
#            assert(nPairs1 == nPairs2)
#            raise Exception
        totalGroundTruth += nPairs1
        if thisFN == 0 and thisFP == 0:
            n_exact += 1
        if q_even:
            N = float(len(refOg)-1)
            totalFN += thisFN/N
            totalFP += thisFP/N
            totalTP += thisTP/N
        else:
            totalFN += thisFN
            totalFP += thisFP
            totalTP += thisTP
        # print("%d\t%d\t%d" % (thisTP, thisFN, thisFP))  
    TP, FP, FN = (totalTP, totalFP, totalFN)
    # print("%d Correct gene pairs" % TP)
    # print("%d False Positives gene pairs" % FP)
    # print("%d False Negatives gene pairs\n" % FN)
    pres = TP/(TP+FP)
    recall = TP/(TP+FN)
    f = 2*pres*recall/(pres+recall)
    print("%0.1f%% F-score" % (100.*f))
    print("%0.1f%% Precision" % (100.*pres))
    print("%0.1f%% Recall\n" % (100.*recall))
    print("%d Orthogroups exactly correct" % n_exact)
    return [100.*f, 100.*pres, 100.*recall]

    
def benchmark(ogs_filename, d_refogs, q_basic_read):
    print("\nReading RefOGs from: %s" % d_refogs)
    ref_ogs = read_refogs(d_refogs)
    ref_ogs_uncertain = read_uncertain_refogs(d_refogs + "low_certainty_assignments/")
    exp_genes = get_expected_genes()
    print("\nReading predicted orthogroups from: %s" % ogs_filename)
    if q_basic_read:
        n_col_skip = get_n_col_skip(ogs_filename)
        pred_ogs = read_orthogroups(ogs_filename, exp_genes, n_col_skip)
    else:
        pred_ogs = read_orthogroups_smart(ogs_filename)
    check_orthogroups(pred_ogs, exp_genes)
    print("\nCalculating benchmarks:")
    # print(os.path.basename(ogs_filename))
    x = calculate_benchmarks_pairwise(ref_ogs, ref_ogs_uncertain, pred_ogs)
    # print("\t".join(map(str, x)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ogs_filename", help="File containing orthogroups.")
    parser.add_argument("-b", "--basic", action="store_true", help="Basic delimited orthogroup file reader")
    args = parser.parse_args()
    d_refogs = "RefOGs" + os.sep
    benchmark(args.ogs_filename, d_refogs, args.basic)
