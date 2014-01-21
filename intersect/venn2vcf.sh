#!/bin/bash

# Converts venn output to the VCF_like, by choosing one of the overlapping
# column data and splitting on a seperator.

## By default the top line contains the file names

usage() {
	echo "Converts venn output to the VCF_like, by splitting on a seperator. It also checks for zygosity consistency"
	echo 
	echo "usage: `basename $0` venn.output --delim=\"S\" [OPTIONS]"
	echo ""
	echo "--delim=N      choose Nth set of data, instead of all"
	echo "--folder=S   set output folder"
	echo ""
	exit -1	
}


if [ $# -lt 1 ] ||  [ $# -gt 4 ]; then
	usage
fi

data_set=-1
sep=-1
ext="_venned.vcf"
folder="./"

mkdir -p $folder

for arg in $*; do
	if [[ "$arg" =~ "--set=" ]]; then
		data_set=`echo $arg | sed 's/--set=//'`
	elif [[ "$arg" =~ "--delim=" ]]; then
		sep=`echo $arg | sed 's/--delim=//'`
	elif [[ "$arg" =~ "--folder=" ]]; then
		folder=`echo $arg | sed 's/--folder=//'`
	fi
done

#Find out how many fields exactly there are
num_fields=`cat $1 | head -5 | awk -F"$sep" '{print NF}'`
names=`cat $1 | head -1 | sed "s/$sep/  /g"`


mkdir -p $folder

#Extract data
if [ "$data_set" != "-1" ];then
	tmp=$RANDOM".data"
	set=$(( $data_set + 1 )) # +1 offset
	awk -F"$sep" -v set=$set '{print $set;}' $1 > $tmp
	name=`echo $names | awk -v ind=$data_set '{print $ind;}'`
	newname=$folder/`echo $name | sed "s/.vcf/$ext/"`
	cat $tmp | egrep -v "^$name" > $newname
	rm $tmp
else
	ind=2
	for name in $names; do
		nnname=`basename $name`
		tmp=$RANDOM".data"
		awk -F"$sep" -v set=$ind '{print $set;}' $1 > $tmp
		newname=$folder/`echo $nnname | sed "s/.vcf/$ext/"`
		cat $tmp | egrep -v "^$name" > $newname
		rm $tmp
		ind=$(( $ind + 1 ))
	done
fi

