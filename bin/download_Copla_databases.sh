#!/usr/bin/env bash

wget https://castillo.dicom.unican.es/zaguan/Copla/Copla_databases_latest.tar
tar -xf Copla_databases_latest.tar

COPLA_DB_DIR=`grep '^COPLA_DB_DIR' copla.ini | cut -f2`

gunzip ${COPLA_DB_DIR}/*.fna.gz

awk -v path=${PWD}/${COPLA_DB_DIR} '{ print path "/" $1 ".fna" }' ${COPLA_DB_DIR}/CoplaDB.lst > ${COPLA_DB_DIR}/CoplaDB.fofn
