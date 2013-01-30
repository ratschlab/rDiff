include bin/rdiff_config.sh

all:	mex

mex:
	echo Entering ./mex
	cd mex ; make octave
	echo Entering ./src/locfit/Source
	cd src/locfit/Source ; make octave

clean:	
	echo Entering ./mex
	cd mex ; make clean
	echo Entering ./examples
	cd examples ; make clean
	echo Entering ./src/locfit/Source
	cd src/locfit/Source ; make clean

example:	
	echo Entering ./examples
	cd examples ; make example

examples:	
	echo Entering ./examples
	cd examples ; make examples
