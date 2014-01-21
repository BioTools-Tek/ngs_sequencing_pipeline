#!/bin/bash

# Extracts and combines multiple gzipped VCF files in a folder
# Batch caller for extractAndCombine_single.sh

usage(){
	echo "`basename $0` <vcf_folder> <output_folder>"
	exit -1
}

if [ $# != 2 ]; then
	usage
fi

vcf_fold=$1
folder=$2

checkStepDone $folder 200

mkdir -p $folder
vcfs=""

for vcf_file in `ls $vcf_fold/*`; do
	vcf=$(extractAndCombine_single.sh $vcf_file)
	printCount $vcf
	mv $vcf $folder/$vcf
	vcf=$folder/$vcf
	vcfs=$vcfs" "$vcf
done

#echo $vcfs
updateLastProcessed $vcfs
