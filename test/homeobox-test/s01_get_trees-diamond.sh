# input variables
n_cpu="8"
alignments="results_trees"
searches="results_searches"

seed_fasta=$1
ref=$2
input_fastas=$3

mkdir -p ${alignments}
mkdir -p ${searches}

# concatenate all proteomes
echo "# Prepare input..."
zcat ${input_fastas}/*.fasta.gz | bioawk -c fastx '{ print ">"$1"\n"$2 }' > ${input_fastas}/all_proteomes.fa
esl-sfetch --index ${input_fastas}/all_proteomes.fa
diamond makedb --in ${input_fastas}/all_proteomes.fa -d ${input_fastas}/all_proteomes.fa --quiet

# DIAMOND search
# alignments
echo "# ${ref} | diamond search"
diamond blastp \
        --more-sensitive \
	--max-target-seqs 100 \
        -d ${input_fastas}/all_proteomes.fa \
        -q ${seed_fasta} \
        -o ${searches}/${ref}.seed.diamond.csv \
        --quiet \
        --threads ${n_cpu}

# extract complete sequences
awk '$11 < 1e-6 { print $2 }' ${searches}/${ref}.seed.diamond.csv | sort -u > ${searches}/${ref}.genes.txt
esl-sfetch -f ${input_fastas}/all_proteomes.fa ${searches}/${ref}.genes.txt > ${alignments}/${ref}.genes.fasta
echo "# ${ref} | found: $(wc -l ${searches}/${ref}.genes.txt)"

# MAFFT, trimming & IQTREE
echo "# ${ref} | Get tree..."
qsub -N tree-${ref} -pe smp ${n_cpu} qsub_alignment.sh ${alignments}/${ref}.genes.fasta ${n_cpu}

