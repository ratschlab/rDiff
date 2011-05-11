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

PROG=`basename $0`
DIR=`dirname $0`

. ${DIR}/../bin/rdiff_config.sh

echo
echo ${PROG}: This program is part of the rDiff version $RDIFF_VERSION.
echo
echo rDiff performs differential expression testing from RNA-Seq measurements.
echo 

if [ -z "$7" -o "$1" == '--help' ];
then
  echo Usage: $0
  echo "   or:" $0 --help
  false
fi 

ANNO_INPUT=${1}
ANNO_FORMAT=${2}
GENES_FN=${3}
BAM_INPUT1=${4}
BAM_INPUT2=${5}
RDIFF_RES_FILE=${6}
TEST_METH=${7}


echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo

if [ "$ANNO_FORMAT" == '0' ]
then
    echo load the genome annotation in GFF3 format and create an annotation object
    if [ ! -f ${GENES_FN} ]
    then
	export PYTHONPATH=$PYTHONPATH:${SCIPY_PATH}
	${PYTHON_PATH} ${DIR}/../tools/ParseGFF.py ${ANNO_INPUT} all all ${GENES_FN}
	${DIR}/../bin/genes_cell2struct ${GENES_FN}
    fi
fi

if [ "$ANNO_FORMAT" == '1' ]
then
    echo load the genome annotation in AGS format
    ln -s ${ANNO_INPUT} ${GENES_FN}
    ${DIR}/../bin/genes_cell2struct ${GENES_FN}
fi

echo
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 2. Differential testing %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo

echo testing genes for differential expression using given alignments
${DIR}/../bin/rdiff ${GENES_FN} ${BAM_INPUT1} ${BAM_INPUT2} ${RDIFF_RES_FILE} ${TEST_METH}

echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo