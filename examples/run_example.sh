
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

echo
echo ${PROG}: This program is part of the rDiff version $RDIFF_VERSION.
echo
echo rDiff performs differential expression testing from RNA-Seq measurements.
echo

if [ -z "$1" -o "$1" == '--help' ];
then
echo Usage: $0 poisson\|param\|nonparam
  echo " or:" $0 --help
  false
fi
if [ "$1" != 'poisson' -a "$1" != 'param' -a "$1" != 'nonparam' ];
then
echo invalid parameter: $1
  false
fi



../bin/rdiff -o results -d data/ -a example_condition_A_replicate_1.bam,example_condition_A_replicate_2.bam -b example_condition_B_replicate_1.bam,example_condition_B_replicate_2.bam -g data/genes_example.gff3 -m $1
