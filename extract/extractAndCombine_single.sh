#!/bin/bash

usage(){
	echo "Extracts VCF archives and combines the indels and snvs into one file">&2
	echo "">&2
	echo "`basename $0` <archive.zip/.tar.gz>">&2
	exit
}

if [ $# != 1 ];then
	usage
fi

file=$1
number=`echo $file | awk -F '/' '{print $NF}' | awk -F'_' '{print $1}'`
number=${number:0:2}
outfile=$number".vcf"

tmp=$RANDOM

#Unpack to tmp
echo -en "\rUnpacking `basename $file | awk -F"---" '{print $1}'`...                ">&2
if [[ "$file" =~ "tar" ]]; then
	direc=`tar xvzf $file | tail -1 2> /dev/null`
	mv $direc $tmp
elif [[ "$file" =~ "zip" ]]; then
	unzip $file -d $tmp > /dev/null 2>&1
fi



#Unpack raw
raw=`ls $tmp/* | grep raw.tgz`
files=`tar xvzf $raw`

echo -en "\rUnpacking `basename $raw`...                      ">&2

headers="headers.txt"
rm $headers $outfile $outfile"_tmp" 2> /dev/null

for f in $files; do
	# annotateType ensures that when we mix indel and snv data,
	# we still know what's what
	if [[ "$f" =~ "indel" ]];then
		annotateType.py $f indel > annotdata
	elif [[ "$f" =~ "snps" ]];then
		annotateType.py $f snv > annotdata
	else
		echo "Problem" >&2
		exit -1
	fi

	egrep "^#" annotdata > header1
	egrep "^#" -v annotdata > data1

	cat header1 >> $headers
	cat data1 >> $outfile"_tmp"

	rm header1 data1 annotdata
done;

rm $files

#Sort Data
echo -en "\rSorting Data...                       ">&2
sort -V $outfile"_tmp" > $outfile"_tmp2"
mv $outfile"_tmp2" $outfile"_tmp"


#Sort headers
echo -en "\rSorting headers...                    ">&2
sort -V $headers | uniq > $headers"_tmp"

egrep "^#CHROM" $headers"_tmp" > chrompos	#Place #CHROM row at bottom, not at top
egrep "^#CHROM" -v $headers"_tmp" > rest

cat rest chrompos > $headers
rm chrompos rest $headers"_tmp"


#Place headers + data
echo -en "\rMaking file.....                      ">&2
cat $headers $outfile"_tmp" > $outfile
rm $headers $outfile"_tmp"
rm -rf $tmp



echo $outfile
