#!/bin/bash

# This script functionally annotates each VCF file, for each gene in the gene list for a single variant:
#
# direction (-|+)
# codon [REF ALT]
# protein [REF ALT]
# mutation [MISSENSE NONSENSE SYNONYMOUS]
# codon change ( c. <num><pos1>_type_<pos2> )
# protein change ( p. <num><REF>_[ALT] )

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

checkStepDone $folder 100 funcannots

files=$(getGroupFiles $group)
mkdir -p $folder

# Functionally annotate
files=$(echo $files | sed 's/ /+/g')
funcannot $files $funcannotopts $folder $folder/"rejects"

updateLastProcessed `ls $folder/*.vcf`

