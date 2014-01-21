#!/bin/bash

usage(){
	echo "`basename $0` <output folder> <group>

 Performs VCF filtering on all files in the group (novel, linkage)
 Batch caller for vcf_filter.py and chrome_comb.py
 
 typical group values are $UNAFF or $AFF
 
 Please initialise_config before running

	"
	exit -1
}

if [ $# -lt 2 ]; then
	usage
fi

folder=$1
group=$2
annotations=$(getGroupFiles $group)

checkStepDone $folder 200 filter

mkdir -p $folder
filters=""

for vcf in $annotations; do
	filt=$(echo `basename $vcf` | sed 's/.vcf/_filter.vcf/')
	filt=$folder/$filt
	
	vcf_filter.py $vcf $vcffopts > $filt
	printCount $filt

	filters=$filters" "$filt
done

#echo $filters
updateLastProcessed $filters

