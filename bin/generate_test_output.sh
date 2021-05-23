#!/usr/bin/env bash

source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate copla

# NZ_CP028167.1
bin/copla.py test/NZ_CP028167.1.fna \
        databases/Copla_RS84/RS84f_sHSBM.pickle \
        databases/Copla_RS84/CoplaDB.fofn \
        test/NZ_CP028167.1.fna_output \
        -a test/NZ_CP028167.1.faa \
        -t circular \
        -k Bacteria \
        -p Proteobacteria \
        -c Gammaproteobacteria \
        -o Enterobacterales \
        -f Enterobacteriaceae \
        -g Escherichia \
        -s 'Escherichia coli' | \
    tee test/NZ_CP028167.1.fna_stdout

tar -zcf test/NZ_CP028167.1.new.tgz -C test NZ_CP028167.1.fna_output/

# NZ_CP028329.1
bin/copla.py test/NZ_CP028329.1.fna \
        databases/Copla_RS84/RS84f_sHSBM.pickle \
        databases/Copla_RS84/CoplaDB.fofn \
        test/NZ_CP028329.1.fna_output \
        -a test/NZ_CP028329.1.faa \
        -t circular \
        -k Bacteria \
        -p Firmicutes \
        -c Bacilli \
        -o Lactobacillales \
        -f Lactobacillaceae \
        -g Lactobacillus \
        -s 'Lactobacillus sp. D1501' | \
    tee test/NZ_CP028329.1.fna_stdout

tar -zcf test/NZ_CP028329.1.new.tgz -C test NZ_CP028329.1.fna_output/
