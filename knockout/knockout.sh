#!/bin/bash

usage(){
	echo "`basename $0` <output folder> <group>"
	exit -1
}

if [ $# != 2 ]; then
	usage
fi

folder=$1
group=$2

files=$(getGroupFiles $group)

checkStepDone $folder 200 knockouts

prefixes=$(split_genes.sh $folder/"lists" $files)
knockouts=$(knockout_prefixes.sh $folder $prefixes)

#echo $knockouts
updateLastProcessed $knockouts
