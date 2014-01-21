#!/bin/bash

# Splits each VCF file into HET, HOM, and CHET components for SNVs and Indels
 

usage(){
	echo "`basename $0` <output folder> <[VCF files]>"
	exit -1
}

if [ $# -lt 2 ]; then
	usage
fi

folder=$1
files="${@:2}"


printline "Splitting HET, HOM, CHET" $txtblu "~" asd
outfiles=""

for file in $files; do
	number=$(echo `basename $file` | awk -F"_" '{print $1}')
	fold=$folder/$number
	outfile=$fold/$number
	
	echo $file >&2
	mkdir -p $fold
	
	hetchethom_split.py $file $outfile $hetchethomopts > $outfile.statistics
	
	mkdir -p $fold/"genes"
	mv $fold/*.genes $fold/"genes"
	
	outfiles=$outfiles" "$outfile
done
echo $outfiles
