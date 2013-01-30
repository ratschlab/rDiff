
#/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
# Copyright (C) 2009-2011 Max Planck Society
#

set -e 

. ../bin/rdiff_config.sh

PROG=`basename $0`

echo
echo ${PROG}: This program is part of the rDiff version $RDIFF_VERSION.
echo
echo rDiff performs differential expression testing from RNA-Seq measurements.
echo 

if [ -z "$1" -o "$1" == '--help' ];
then
  echo Usage: $0 poisson\|param\|nonparam
  echo "   or:" $0 --help
  false
fi
if [ "$1" != 'poisson' -a "$1" != 'param' -a "$1" != 'nonparam' ];
then
  echo invalid parameter: $1
  false
fi

RDIFF_INPUT_DIR=data/
#FASTA_INPUT=data/TAIR10.fasta
GFF3_INPUT="${RDIFF_INPUT_DIR}/genes_example.gff3"
#SAM_INPUT1=data/c_elegans_WS200-I-regions-SRX001872.sam
#SAM_INPUT2=data/c_elegans_WS200-I-regions-SRX001875.sam
                                                                                                                                                                                                                                                                                                                                                              
BAM_INPUT1="example_condition_A_replicate_1.bam,example_condition_A_replicate_2.bam"
BAM_INPUT2="example_condition_B_replicate_1.bam,example_condition_B_replicate_2.bam"


EXP=Arti
TEST_METH=$1
if [ "$1" == 'poisson' ];
then
    echo Note: Running this script takes about 5 minute \(on a single CPU\).
    TEST_METH_NAME=poisson
fi
if [ "$1" == 'param' ];
then
    echo Note: Running this script takes about 5 minute \(on a single CPU\).
    TEST_METH_NAME=parametric
fi
if [ "$1" == 'nonparam' ];
then
    echo Note: Running this script takes about 60 minutes \(on a single CPU\).
    TEST_METH_NAME=nonparametric
fi

RESULTDIR=./results-$1
mkdir -p $RESULTDIR

echo All results can be found in $RESULTDIR
echo

echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo

GENES_RES_DIR=$RESULTDIR/genes.anno
mkdir -p $GENES_RES_DIR
GENES_FN=${GENES_RES_DIR}/genes.mat

echo 1a. load the genome annotation in GFF3 format, create an annotation object #\[Log file in ${GENES_RES_DIR}/elegans-gff2anno.log\]
if [ ! -f ${GENES_FN} ]
then
    export PYTHONPATH=$PYTHONPATH:${SCIPY_PATH}
echo "${RDIFF_PYTHON_PATH} -W ignore::FutureWarning ../tools/ParseGFF.py ${GFF3_INPUT} ${GENES_FN} "
    ${RDIFF_PYTHON_PATH} -W ignore::FutureWarning ../tools/ParseGFF.py ${GFF3_INPUT} ${GENES_FN} #> ${RESULTDIR}/elegans-gff2anno.log
    ../bin/genes_cell2struct ${GENES_FN}
fi


#echo 1b. convert the alignments in SAM format to BAM format
#../tools/./sam_to_bam.sh ${SAMTOOLS_DIR} data ${FASTA_INPUT} ${SAM_INPUT1}
#../tools/./sam_to_bam.sh ${SAMTOOLS_DIR} data ${FASTA_INPUT} ${SAM_INPUT2}

echo
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 2. Differential testing %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo

RDIFF_RES_DIR=$RESULTDIR
mkdir -p $RDIFF_RES_DIR
RDIFF_RES_FILE=$RDIFF_RES_DIR/P_values_rDiff_${TEST_METH_NAME}.tab

echo "../bin/rDiff ${GENES_FN} ${BAM_INPUT1} ${BAM_INPUT2} ${RDIFF_RES_FILE} ${TEST_METH} ${RESULTDIR} ${RDIFF_INPUT_DIR} > ${RESULTDIR}/example-rdiff.log"
echo testing genes for differential expression using given alignments \(log file in ${RESULTDIR}/example-rdiff.log\)
../bin/rdiff ${GENES_FN} ${BAM_INPUT1} ${BAM_INPUT2} ${RDIFF_RES_FILE} ${TEST_METH} ${RESULTDIR} ${RDIFF_INPUT_DIR} #> ${RESULTDIR}/example-rdiff.log

echo
echo Testing result can be found in $RDIFF_RES_FILE
echo

echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo