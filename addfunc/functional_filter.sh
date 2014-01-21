#!/bin/bash

# This script filters pre-functionally annotated VCF files (functional.sh)
# against keywords (e.g omit SYNONYMOUS) and other parameters (e.g
# omit isoforms and splice sites).

# Then it is used in conjunction with common.sh script to get common
# variants for the remaining genes.

# Then later used with genetables.sh script to display gene information
# per file.


usage(){
	echo "`basename $0` <output_folder> <group>"
	echo ""
	echo "Filter keywords are set by editing the keywops variable in load_config.sh" 
	echo "Isoform and Splice removal (defaults) are set by editing the isointrosplops variable in the same file"
	echo ""
	echo "Please run . initialise_config.sh first to set variables."
	exit
}

[ $# != 2 ] && usage

folder=$1
group=$2

files=$(getGroupFiles $group)

isodir=$folder/"1_isoform_intron_splice_out"
fildir=$folder/"2_keyword_filter_in"

mkdir -p $isodir # $fildir

outputfiles=""

echo "FILES:"$files >&2

for file in $files; do
	num=$(getNum $file);
	
	outfile1=$isodir/$num"_isointospl.vcf"
	outfile2=$fildir/$num"_funcannot_filter.vcf"
	
	outfile=""
	
	#Filter out splice or introns first (keep only exons)
	if [ "$isointrosplops" != "" ];then
		mkdir -p $isodir
		filter_funcannot_isointrospl.py $file $allopts $isointrosplops > $outfile1
		file=$outfile1			# to be read by keywops
		outfile=$outfile1		# last outted file
	fi
	
# Deprecated by gene_tables_HTML.py which enables interactive filtering on Keyword (well... synonymous anyway).
	# Filter on keyword first
#	if [ "$keywops" != "" ]; then
#		mkdir -p $fildir
#		filter_funcannot.py $file $allopts $keywops > $outfile2
#		outfile=$outfile2		# last outted file
#	fi
	outputfiles=$outputfiles" "$outfile
done

updateLastProcessed $outputfiles


#Now get remainder variants from remaining overlapping genes
commonfold=$folder/"3_overlapping_common"
$(common.sh $commonfold $group)


