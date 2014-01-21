#!/bin/bash

if [ "$delim" = "" ];then
	echo "run . initialise_config.sh first"
	exit
fi

usage(){
	echo "`basename $0` <output fold> <[affected_VCFs]> <[unaffected_VCFs]>
	
	where individual files are delimited by $delim for each seperate file list
	"
	exit
}

if [ $# != 3 ];then
	usage
fi

folder=$1
affecteds=`echo $2 | sed s"/$delim/ /g"`
unaffecteds=`echo $3 | sed s"/$delim/ /g"`

genefold=$folder/"a_affected_genelists"
matchfold=$folder/"b_matched_unaffecteds_with_all_affected_genes"

mkdir -p $genefold
mkdir -p $matchfold

#Produce an all encompassing gene list from all affecteds
gene_files=""
for aff in $affecteds; do
	genes=$genefold/$(echo `basename $aff` | sed s'/.vcf/.genes/')
	extract_genes_vcf.py $geneidopts $aff > $genes	# Do NOT include isoforms or exoforms -- need to be as general as possible
	gene_files=$gene_files" "$genes
done

#Sort and uniq all genelists into one combined
combined=$folder/"all_affected_genes.list"
cat $gene_files | sort | uniq > $combined

# Matched the unaffected to those genes
matched=""
for unaff in $unaffecteds; do
	name=$matchfold/$(echo `basename $unaff` | awk -F "_" '{print $1}')"_matched_genes.vcf"
	match_gene_lines.py $unaff $combined $matchglopts > $name   ## no exons or isoforms, but must give geneID
	matched=$matched" "$name
done

echo $matched
