#!/bin/bash

# Search for pedfile
export pedfile="pedfile.pro"
if [ ! -f $pedfile ];then
	echo -e "\nPlease include a pedfile.pro file in the directory\n" >&2
	exit -1
fi

# Loads ext modules
. set_path.sh ~/bin/   	#.bashrc
. load_config.sh
. load_folders.sh
. load_functions.sh
