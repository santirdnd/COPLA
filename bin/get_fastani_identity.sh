#!/usr/bin/env bash

if [ $# -ne 2 ]; then
    echo -e 'Usage: get_ani_identity.sh QUERY.fna.gz REFLIST.lst'
    exit 2
fi

source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate /fernando/envs/fastani

THREADS=48
MINFRAC=0.5
FRAGLEN=1500
KMERLEN=16

QRYFILE=$1
REFLIST=$2
OUTFILE="$1.fastani.tsv"
LOGFILE="$1.fastani.log"

fastANI -q ${QRYFILE} --rl ${REFLIST} -o ${OUTFILE} --minFraction ${MINFRAC} -t ${THREADS} -k ${KMERLEN} --fragLen ${FRAGLEN} 1>${LOGFILE} 2>&1
