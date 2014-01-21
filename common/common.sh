#!/bin/bash
# Keep variants in all files that exist within the common genelist
# across all files

usage(){
	echo "`basename $0` <output folder> <group>
	"
	exit
}

if [ $# != 2 ];then
	usage
fi

folder=$1
group=$2

files=$(getGroupFiles $group)
mkdir -p $folder

goptname="opts_"`echo "$commonopts" | sed 's/--//g' | sed 's/ /_/g'`
printline " Options: $gopts " $txtpur "." asd

g_fold=$folder/$goptname/"genes"
m_fold=$folder/$goptname/"match"

mkdir -p $g_fold $m_fold;

glist=""
for file in $files; do
	gene=$g_fold/$(echo `basename $file` | sed 's/.vcf/.genes/')
	extract_genes_vcf.py $geneidopts $file $commonopts > $gene
	glist=$glist" "$gene
done
common=$folder/$goptname/"common_genelist_$goptname.txt"
overlap_genes.py $glist > $common

#Grab matching
newnames=""
for file in $files; do
	newname=$m_fold/$(getNum $file)"_matched_common.vcf"
	match_gene_lines.py $file $common $matchglopts > $newname
	newnames=$newnames" "$newname
done

updateLastProcessed $newnames
