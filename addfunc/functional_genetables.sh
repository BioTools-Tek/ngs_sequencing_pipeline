#!/bin/bash

# This is a wrapper script for genetables.sh which requires functional_filter.sh to run first.

usage(){
	echo "`basename $0` <output_folder> <group>"
	echo ""
	echo "Please run . initialise_config.sh first"
	exit
}

if [ $# -lt 2 ]; then
   usage
fi

folder=$1
group=$2

#checkStepDone $folder

# Apply filtering
$(functional_filter.sh $folder $group)

# Produce gene tables
newfold=$folder/"4_genetables"
echo $(genetables.sh $newfold $(getGroupFiles $group))  # echo's genetable.csv

