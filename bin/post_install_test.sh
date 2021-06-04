#!/usr/bin/env bash

main(){
    check_bash
    check_conda
    check_perl
    check_python3

    check_graph_tool
    check_numpy
    check_pandas

    check_blastn
    check_prodigal
    check_hmmscan
    check_plasmidfinder
    check_macsyfinder

    check_ruby
    check_parallel
    check_ani_rb

    test_plasmid_NZ_CP028167
    test_plasmid_NZ_CP028329
}

version_check(){
    # Based on https://unix.stackexchange.com/questions/285924/how-to-compare-a-programs-version-in-a-shell-script
    VERSION=$1
    REQUIRED=$2

    if [ "$(printf '%s\n' "${VERSION}" "${REQUIRED}" | sort -V | head -n1)" = "${REQUIRED}" ]; then
        echo "  OK - Current version: ${VERSION}"
        return 0
    else
        echo '  Warning! Please verify the installed version'
        echo "    Current version: ${VERSION}"
        echo "    Minimum version required: ${REQUIRED}"
        return 1
    fi
}

check_bash(){
    EXE='bash'
    echo "Checking for ${EXE} ..."

    EXE_PATH=`which ${EXE}`
    if [ -z "${EXE_PATH}" ]; then
        echo "  Error! ${EXE} is not in your PATH"
    else
        VERSION=`bash --version | head -n1 | cut -d' ' -f4`
        REQUIRED='4.0'
        version_check ${VERSION} ${REQUIRED}
    fi
}

check_conda(){
    EXE='conda'
    echo "Checking for ${EXE} ..."

    EXE_PATH=`which ${EXE}`
    if [ -z "${EXE_PATH}" ]; then
        echo "  Error! ${EXE} is not in your PATH"
    else
        VERSION=`conda --version | cut -d' ' -f2`
        REQUIRED='4.5'
        version_check ${VERSION} ${REQUIRED}
    fi

    CONDA_SHELL_INT="$(dirname "${CONDA_EXE%/*}")"/etc/profile.d/conda.sh
    if [ ! -e "${CONDA_SHELL_INT}" ]; then
        echo "  Error! Please check how conda integrates with your shell"
        echo "         For copla to activate the macsyfinder environment it is assumed that ${CONDA_SHELL_INT} can be sourced"
    fi
}

check_perl(){
    EXE='perl'
    echo "Checking for ${EXE} ..."

    EXE_PATH=`which ${EXE}`
    if [ -z "${EXE_PATH}" ]; then
        echo "  Error! ${EXE} is not in your PATH"
    else
        VERSION=`perl -v | head -n2 | tail -n+2 | cut -d' ' -f9 | tr -d [\(\)]`
        REQUIRED='5.22'
        version_check ${VERSION} ${REQUIRED}
    fi
}

check_python3(){
    EXE='python3'
    echo "Checking for ${EXE} ..."

    EXE_PATH=`which ${EXE}`
    if [ -z "${EXE_PATH}" ]; then
        echo "  Error! ${EXE} is not in your PATH"
    else
        VERSION=`python3 -V | cut -d' ' -f2`
        REQUIRED='3.8'
        version_check ${VERSION} ${REQUIRED}
    fi
}

check_graph_tool(){
    PYTHON_MODULE='graph_tool'
    echo "Checking for ${PYTHON_MODULE} ..."

    MODULE_LOAD=`python3 -c "import ${PYTHON_MODULE}" 2>&1`
    if [ -z "${MODULE_LOAD}" ]; then
        VERSION=`python3 -c "import ${PYTHON_MODULE}; print(${PYTHON_MODULE}.__version__)" | cut -d' ' -f1`
        REQUIRED='2.33'
        version_check ${VERSION} ${REQUIRED}
    else
        echo "  Error! ${PYTHON_MODULE} is not in your PATH"
    fi
}

check_numpy(){
    PYTHON_MODULE='numpy'
    echo "Checking for ${PYTHON_MODULE} ..."

    MODULE_LOAD=`python3 -c "import ${PYTHON_MODULE}" 2>&1 `
    if [ -z "${MODULE_LOAD}" ]; then
        VERSION=`python3 -c "import ${PYTHON_MODULE}; print(${PYTHON_MODULE}.__version__)"`
        REQUIRED='1.19.1'
        version_check ${VERSION} ${REQUIRED}
    else
        echo "  Error! ${PYTHON_MODULE} is not in your PATH"
    fi
}

check_pandas(){
    PYTHON_MODULE='pandas'
    echo "Checking for ${PYTHON_MODULE} ..."

    MODULE_LOAD=`python3 -c "import ${PYTHON_MODULE}" 2>&1`
    if [ -z "${MODULE_LOAD}" ]; then
        VERSION=`python3 -c "import ${PYTHON_MODULE}; print(${PYTHON_MODULE}.__version__)"`
        REQUIRED='1.1.0'
        version_check ${VERSION} ${REQUIRED}
    else
        echo "  Error! ${PYTHON_MODULE} is not in your PATH"
    fi
}

check_blastn(){
    EXE='blastn'
    echo "Checking for ${EXE} ..."

    EXE_PATH=`which ${EXE}`
    if [ -z "${EXE_PATH}" ]; then
        echo "  Error! ${EXE} is not in your PATH"
    else
        VERSION=`blastn -version | head -n1 | cut -d' ' -f2`
        REQUIRED='2.9.0'
        version_check ${VERSION} ${REQUIRED}
    fi
}

check_prodigal(){
    EXE='prodigal'
    echo "Checking for ${EXE} ..."

    EXE_PATH=`which ${EXE}`
    if [ -z "${EXE_PATH}" ]; then
        echo "  Error! ${EXE} is not in your PATH"
    else
        VERSION=`prodigal -v 2>&1 | head -n2 | tail -n+2 | cut -d' ' -f2 | sed 's/^V//' | sed 's/:$//'`
        REQUIRED='2.6.3'
        version_check ${VERSION} ${REQUIRED}
    fi
}

check_hmmscan(){
    EXE='hmmscan'
    echo "Checking for ${EXE} ..."

    EXE_PATH=`which ${EXE}`
    if [ -z "${EXE_PATH}" ]; then
        echo "  Error! ${EXE} is not in your PATH"
    else
        VERSION=`hmmscan -h | head -n2 | tail -n+2 | cut -d' ' -f3`
        REQUIRED='3.1'
        version_check ${VERSION} ${REQUIRED}
    fi
}

check_plasmidfinder(){
    EXE='plasmidfinder.py'
    echo "Checking for ${EXE} ..."

    EXE_PATH=`which ${EXE}`
    if [ -z "${EXE_PATH}" ]; then
        echo "  Error! ${EXE} is not in your PATH"
    else
        echo '  OK'
    fi
}

check_macsyfinder(){
    EXE='macsyfinder'
    echo "Checking for ${EXE} ..."

    EXE_CONDA=conda
    EXE_CONDA_PATH=`which ${EXE_CONDA}`
    if [ -z "${EXE_CONDA_PATH}" ]; then
        EXE_PATH=`which ${EXE}`
        if [ -z "${EXE_PATH}" ]; then
            echo "  Error! ${EXE} is not in your PATH"
        else
            VERSION=`macsyfinder --version 2>&1 | head -n1 | cut -d' ' -f2`
            REQUIRED='1.0.5'
            version_check ${VERSION} ${REQUIRED}
        fi
    else
        source "$(dirname "${CONDA_EXE%/*}")"/etc/profile.d/conda.sh
        conda activate macsyfinder

        EXE_PATH=`which ${EXE}`
        if [ -z "${EXE_PATH}" ]; then
            echo "  Error! ${EXE} is not in your PATH. Please verify the macsyfinder conda environment"
        else
            VERSION=`macsyfinder --version 2>&1 | head -n1 | cut -d' ' -f2`
            REQUIRED='1.0.5'
            version_check ${VERSION} ${REQUIRED}
        fi

        conda deactivate
        conda activate copla
    fi
}

check_ruby(){
    EXE='ruby'
    echo "Checking for ${EXE} ..."

    EXE_PATH=`which ${EXE}`
    if [ -z "${EXE_PATH}" ]; then
        echo "  Error! ${EXE} is not in your PATH"
    else
        VERSION=`ruby -v | cut -d' ' -f2`
        REQUIRED='2.3'
        version_check ${VERSION} ${REQUIRED}
    fi
}

check_parallel(){
    EXE='parallel'
    echo "Checking for ${EXE} ..."

    EXE_PATH=`which ${EXE}`
    if [ -z "${EXE_PATH}" ]; then
        echo "  Error! ${EXE} is not in your PATH"
    else
        VERSION=`parallel --version | head -n1 | cut -d' ' -f3`
        REQUIRED='20161222'
        version_check ${VERSION} ${REQUIRED}
    fi
}

check_ani_rb(){
    EXE='ani.rb'
    echo "Checking for ${EXE} ..."

    EXE_PATH='bin/${EXE}'
    if [ -z "${EXE_PATH}" ]; then
        echo "  Error! ${EXE} is not in your PATH"
    else
        echo '  OK'
    fi
}

test_plasmid_NZ_CP028167(){
    SEQ_FNA='test/NZ_CP028167.1.fna'
    SEQ_FAA='test/NZ_CP028167.1.faa'

    COPLA_DB=`grep '^COPLA_DB_DIR' copla.ini | cut -f2`
    OUTPUT_DIR=${SEQ_FNA}_output
    OUTPUT=${SEQ_FNA}_stdout
    TMP_FILE=${SEQ_FNA}_tmp

    echo 'Checking COPLA with plasmid NZ_CP028167.1 ...'
    echo '=========='
    bin/copla.py ${SEQ_FNA} ${COPLA_DB}/RS84f_sHSBM.pickle ${COPLA_DB}/CoplaDB.fofn ${OUTPUT_DIR} \
        -a ${SEQ_FAA} -t circular -k Bacteria -p Proteobacteria -c Gammaproteobacteria -o Enterobacterales \
        -f Enterobacteriaceae -g Escherichia -s 'Escherichia coli' | tee ${TMP_FILE}
    echo '=========='
    DIFF_OUT=`diff ${OUTPUT} ${TMP_FILE}`
    if [ -z "${DIFF_OUT}" ]; then
        echo '  OK - Output for plasmid NZ_CP028167.1 is as expected'
    else
        echo '  Error! Output for plasmid NZ_CP028167.1 differs. Please verify COPLA installation'
    fi
    rm ${TMP_FILE}
}

test_plasmid_NZ_CP028329(){
    SEQ_FNA='test/NZ_CP028329.1.fna'

    COPLA_DB=`grep '^COPLA_DB_DIR' copla.ini | cut -f2`
    OUTPUT_DIR=${SEQ_FNA}_output
    OUTPUT=${SEQ_FNA}_stdout
    TMP_FILE=${SEQ_FNA}_tmp

    echo 'Checking COPLA with plasmid NZ_CP028329.1 ...'
    echo '=========='
    bin/copla.py ${SEQ_FNA} ${COPLA_DB}/RS84f_sHSBM.pickle ${COPLA_DB}/CoplaDB.fofn ${OUTPUT_DIR} \
        -t circular -k Bacteria -p Firmicutes -c Bacilli -o Lactobacillales -f Lactobacillaceae \
        -g Lactobacillus -s 'Lactobacillus sp. D1501' | tee ${TMP_FILE}
    echo '=========='
    DIFF_OUT=`diff ${OUTPUT} ${TMP_FILE}`
    if [ -z "${DIFF_OUT}" ]; then
        echo '  OK - Output for plasmid NZ_CP028329.1 is as expected'
    else
        echo '  Error! Output for plasmid NZ_CP028329.1 differs. Please verify COPLA installation'
    fi
    rm ${TMP_FILE}
}

main
