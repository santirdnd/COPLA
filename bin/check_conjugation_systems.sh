#!/usr/bin/env bash

if [ $# -ne 5 ]; then
    echo 'Usage: check_conjugation_systems.sh QUERY.faa OUTPUT.path TYPE TOPOLOGY DATABASE.path'
    echo '    TYPE: ordered_replicon or unordered_replicon'
    echo '    TOPOLOGY: circular or linear'
    exit 2
fi

source "$(dirname "${CONDA_EXE%/*}")"/etc/profile.d/conda.sh
conda activate macsyfinder

QRYFILE=$1
OUT_DIR=$2
SEQTYPE=$3
SEQTOPO=$4
DB_PATH=$5
DB_DEFS=${DB_PATH}/definitions/
DB_PRFS=${DB_PATH}/profiles/

THREADS=20
OUTFILE=${OUT_DIR}/'results_tab.tsv'
OUTFILE_REPORT=${OUT_DIR}/'results_tab.report.tsv'
OUTFILE_SUMMARY=${OUT_DIR}/'results_tab.summary.tsv'

mkdir -p ${OUT_DIR}
for conj_type in typeF typeB typeC typeFATA typeFA typeG typeI typeT; do
    macsyfinder ${conj_type} \
        -w ${THREADS} \
        -d ${DB_DEFS} \
        -p ${DB_PRFS} \
        --sequence-db ${QRYFILE} \
        --db-type ${SEQTYPE} \
        --replicon-topology ${SEQTOPO} \
        -o ${OUT_DIR}/${conj_type}
done

# Plasmids with no MPF will have no result files. Must check if there is at least a file
if compgen -G "${OUT_DIR}/*/macsyfinder.tab" > /dev/null; then
    awk 'ORS=NR%2?"\t":"\n"' ${OUT_DIR}/*/macsyfinder.tab | cut -f2- > ${OUTFILE}
fi
if compgen -G "${OUT_DIR}/*/macsyfinder.report" > /dev/null; then
    awk 'NR==1||FNR>1' ${OUT_DIR}/*/macsyfinder.report > ${OUTFILE_REPORT}
fi
if compgen -G "${OUT_DIR}/*/macsyfinder.summary" > /dev/null; then
    awk 'NR==1||FNR>1' ${OUT_DIR}/*/macsyfinder.summary > ${OUTFILE_SUMMARY}
fi

rm ${QRYFILE}.* $(dirname ${QRYFILE})/formatdb.err
