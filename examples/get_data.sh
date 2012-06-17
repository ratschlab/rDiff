#!/bin/bash

set -e

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % rDiff examples: get_data.sh %
echo % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo 
echo This script gets reads and the genome sequence 
echo for one C. elegans example.
echo 

export DATA_DIR=data

if [ -d $DATA_DIR ]
then
	echo Data directory ./$DATA_DIR already exists
	echo \(remove it to run this script again\)
	exit -1
fi

echo Downloading rDiff example data from FTP server ...
wget -c ftp://ftp.tuebingen.mpg.de/fml/bohnert/rDiff/rdiff_examples.tar.bz2
echo uncompressing ...
tar -xjf rdiff_examples.tar.bz2
mv rdiff_examples $DATA_DIR
echo
echo -n Done.
echo
