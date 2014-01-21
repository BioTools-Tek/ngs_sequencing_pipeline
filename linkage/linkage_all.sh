#!/bin/bash

usage(){
	echo "Filters based upon a bed file (i.e. linkage)

	`basename $0` <output folder> <linkage.bed> <group>
or	`basename $0` <output folder> <linkage.bed> <[FILES]>
" >&2
	exit
}

folder=$1
bed=$2
vcfs=""

if [ $# = 3 ];then
	group=$3
	checkStepDone $folder 100 linkage  # Has to look 4 subdirs to find last updated files
	vcfs=$(getGroupFiles $group)

elif [ $# -gt 3 ];then
	vcfs="${@:3}"
else
	usage
fi


[ "$bed" = "no" ] || [ "$bed" = "" ] && exit

mkdir -p $folder

for vcf in $vcfs; do
	id=$(getNum $vcf)
	newname=$folder/$id"_linkage.vcf"
	
	chrome_comb.py $vcf $bed > $newname
done

updateLastProcessed `ls $folder/*.vcf`   #needed by common.sh to know what files to check
common.sh $folder/"genes" $group

