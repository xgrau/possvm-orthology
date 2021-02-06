# input variables
n_cpu="8"
hmm_folder="hmm_profiles_weak/"
input_fastas="proteomes/"
alignments="hmm_trees"
searches="hmm_searches"

mkdir -p ${alignments}
mkdir -p ${searches}

# functions
function do_hmmsearch {

    # input
    local hmm=$1
    local out=$2
    local fas=$3
    local cpu=$4
    local thr=$5

    if [ $thr == "GA" ] ; then
    # search with gathering threshold
    hmmsearch \
        --domtblout ${out}.domtable \
        --cut_ga \
        --cpu ${cpu} \
        ${hmm} \
        ${fas} 1> /dev/null
    else
    # search with domain evalue
    hmmsearch \
        --domtblout ${out}.domtable \
        --domE ${thr} \
        --cpu ${cpu} \
        ${hmm} \
        ${fas} 1> /dev/null
    fi

    # clean domtable
    grep -v "^#" ${out}.domtable \
    | awk 'BEGIN { OFS="\t" } { print $1,$18,$19,$4,$5,$12 }' \
    > ${out}.domtable.csv

}

function do_alitrimphy {

	local i=$1
	local c=$2
	local il=$3
	local ilt=$4
	local tre=$5

	# align & trim
	if [ -s $il ] ; then
	echo "$il already exists, go to next step"
	else
	mafft --genafpair --thread $c --reorder --maxiterate 1000 $i > $il 2> /dev/null
	fi
	if [ -s $ilt ] ; then
        echo "$ilt already exists, go to next step"
	else
	clipkit $il -m kpic-gappy -o $ilt -g 0.7 2> /dev/null 1> /dev/null
	fi

	# iqtree
	iqtree \
		-s $ilt \
		-m TEST \
		-mset LG,WAG,JTT \
		-nt AUTO \
		-ntmax $c \
		-bb 1000 \
		-pre $tre \
		-cptime 1800 2> /dev/null 1> /dev/null

}

# concatenate all proteomes
echo "# Prepare input..."
zcat ${input_fastas}/*.fasta.gz > ${input_fastas}/all_proteomes.fa
esl-sfetch --index ${input_fastas}/all_proteomes.fa

if [ -z $1 ]; then
	hmm_list=${hmm_folder}/*.hmm
else
	hmm_list=$1
fi

# do HMM searches for each refOG seed HMM
for hmm in ${hmm_list} ; do

	# HMM name
	ref=$(basename ${hmm%%.hmm})

	# HMM search
	echo "# ${ref} | HMM search"
	do_hmmsearch ${hmm} ${searches}/${ref}.hmmsearch ${input_fastas}/all_proteomes.fa ${n_cpu} 1e-15
	# extract complete sequences
	cut -f1 ${searches}/${ref}.hmmsearch.domtable.csv | sort -u > ${searches}/${ref}.hmmsearch.genes.txt
	esl-sfetch -f ${input_fastas}/all_proteomes.fa ${searches}/${ref}.hmmsearch.genes.txt > ${alignments}/${ref}.genes.fasta
	echo "# ${ref} | found: $(wc -l ${searches}/${ref}.hmmsearch.genes.txt)"
	echo "# ${ref} | refOG: $(fgrep -w ${ref%%.*} refOGs.csv | cut -f2 | sort -u | wc -l)"
	echo "# ${ref} | refOG in found: $(comm -12 ${searches}/${ref}.hmmsearch.genes.txt <(fgrep -w ${ref%%.*} refOGs.csv | cut -f2 | sort -u) | wc -l)"
	if [ $(fgrep -w ${ref%%.*} refOGs.csv | cut -f2 | sort -u | wc -l) -ne $(comm -12 ${searches}/${ref}.hmmsearch.genes.txt <(fgrep -w ${ref%%.*} refOGs.csv | cut -f2 | sort -u) | wc -l) ] ; then
	echo "# ${ref} | missing: $(comm -13 ${searches}/${ref}.hmmsearch.genes.txt <(fgrep -w ${ref%%.*} refOGs.csv | cut -f2 | sort -u) | tr '\n' ',')"
	fi

done


# do MSA and trees for each HMM
for hmm in ${hmm_list} ; do

	# HMM name
	ref=$(basename ${hmm%%.hmm})

	# MAFFT, trimming & IQTREE
	echo "# ${ref} | Get tree..."
	do_alitrimphy ${alignments}/${ref}.genes.fasta ${n_cpu} ${alignments}/${ref}.genes.l.fasta ${alignments}/${ref}.genes.lt.fasta ${alignments}/${ref}.genes.iqtree

done

