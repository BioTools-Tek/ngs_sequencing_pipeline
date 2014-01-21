#!/bin/bash

# Combines the HET HOM CHET gene names to form knockout and 
# non knockout gene lists
 

usage(){
	echo "`basename $0` <output folder> <[list of prefixes]>
	
where a prefix denotes where to find the HET HOM CHET files for an individual VCF file"
	exit -1
}

if [ $# -lt 2 ]; then
	usage
fi

folder=$1
prefixes="${@:2}"

mkdir -p $folder


printline " Combining KO genes " $txtblu "~" asd
knockouts=""

for pref in $prefixes; do
	number=$(echo `basename $pref` | awk -F"_" '{print $1}')
	files=`ls $pref*.vcf`

	headers=$folder/$number"_headers"
	data=$folder/$number"_data"
	outfile=$folder/$number"_knockouts.vcf"
	rm $headers $data 2>/dev/null
	
	for f in $files; do
		if [[ $f =~ \.het\.vcf ]];then
			continue
		fi

		egrep "^#" $f >> $headers
		egrep "^#" -v $f >> $data

	done


	## Sort headers
	sort -V $headers | uniq > $headers"_tmp"
	egrep "^#CHROM" -v $headers"_tmp" > $headers"_rest"	#Place #CHROM row at bottom, not at top
	egrep "^#CHROM" $headers"_tmp" > $headers"_chrompos"

	cat $headers"_rest" $headers"_chrompos" > $headers
	rm $headers"_rest" $headers"_chrompos" $headers"_tmp"

	## Sort and uniq data
	sort -V $data | uniq > $data"_tmp"
	mv $data"_tmp" $data

	## Combine
	cat $headers $data > $outfile
	rm $headers $data
	
	printCount $outfile
	
	knockouts=$knockouts" "$outfile

done

echo $knockouts
