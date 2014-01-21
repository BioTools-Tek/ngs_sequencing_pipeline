#!/bin/bash

usage(){
	echo -e "
   `basename $0` <output_fold> <group>

This script takes in a list of pre-annotated pre-filtered files (removing introns) (i.e. run functional_filter.sh first before thi script) and displays them in a table format showing how many mutations occur per gene per patient, with gene on the rows and patients as columns

Please edit the genetableflags variable in load_config.sh script to specify filtering by genename or exon, include frequencies, and to print types of mutations rather than nomenclature.
Then run . initialise_config.sh before running this script.
"
	exit
}

[ $# != 2 ] && usage

folder=$1
group=$2
files=$(getGroupFiles $group)

# natural sort files
files=$(echo $files | sed 's/ /\n/g' | sort -g)

mkdir -p $folder 
table=$folder/"gene_table.html"

gene_tables_HTML.py $genetableopts $genetable_flags $files > $table

nohup firefox $table &
