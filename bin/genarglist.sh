#/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
# Copyright (C) 2009-2010 Max Planck Society
#

until [ -z $1 ] ; do
	if [ $# != 1 ];
	then
		echo -n "$1:"
	else
		echo -n "$1"
	fi
	shift
done
