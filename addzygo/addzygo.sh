#!/bin/bash

usage(){
	echo "`basename $0` <output folder> <group>
	
Appends zygosity to all files
Batch caller for bamzygo
"
	exit -1
}

if [ $# != 2 ]; then
	usage
fi

folder=$1
group=$2

files=$(getGroupFiles $group)

#checkStepDone $folder 200 zygosity

zygosity_files=""
mkdir -p $folder

for vcf in $files; do
	number=`echo $vcf | awk -F '/' '{print $NF}' | awk -F'_' '{print $1}'`
	number=${number:0:2}

	unknowns=$folder/$(echo `basename $vcf` | sed 's/.vcf/_nobamdata.vcf/' )
	zygoes=$folder/$(echo `basename $vcf` | sed 's/.vcf/_zygosity.vcf/' )

	if [ -f $zygoes ]; then
		printCount $zygoes
	else
		echo $zygoes >&2
		bamzygo "$bam_fold/$number""_"*.bam $vcf $unknowns $bamzyopts > $zygoes
	fi
	zygosity_files=$zygosity_files" "$zygoes
done

#echo $zygosity_files
updateLastProcessed $zygosity_files

