include bin/rdiff_config.sh

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2013 Gunnar Raetsch, Regina Bohnert, Philipp Drewe
# Copyright (C) 2009-2013 Max Planck Society, Sloan-Kettering Institute
#

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
