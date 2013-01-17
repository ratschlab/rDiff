#!/bin/bash

set -e

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % rDiff examples: get_data.sh %
echo % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo 
echo This script gets reads and the genome sequence 
echo for one A.thaliana example.
echo 

export DATA_DIR=data

if [ -d $DATA_DIR ]
then
	echo Data directory ./$DATA_DIR already exists
	echo \(remove it to run this script again\)
	exit -1
fi

echo Downloading rDiff example data from FTP server ...
wget -c http://cbio.mskcc.org/public/raetschlab/user/drewe/rdiff/example_data.tar.gz
echo uncompressing ...
tar -xzf example_data.tar.gz
mv example_data $DATA_DIR
echo
echo -n Done.
echo
