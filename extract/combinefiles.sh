#!/bin/bash

if [ $# -lt 3 ]; then
	echo "Merges unique data from two or more VCF file

	usage: `basename $0` <newname.vcf> <[Files.VCF]>" >&2
	exit
fi

files=${@:2}

outfile=$1
headers=$outfile"._headers"
chrompos=$outfile"_chrompos"
rest=$outfile"_rest"
offs="\t\t\t\t\t\t\t\t"

echo -en "\rCombining:" >&2
for f in $files; do
	echo -en " `basename $f`" >&2

	header1=$f"._headers"
	data1=$f"._data"

	egrep "^#" $f > $header1
	egrep "^#" -v $f > $data1

	cat $header1 >> $headers
	cat $data1 >> $outfile"_tmp"

	rm $header1 $data1
done;

#Sort Data
echo -en "\r$offs Data ">&2
sort -V $outfile"_tmp" > $outfile"_tmp2"
mv $outfile"_tmp2" $outfile"_tmp"


#Sort headers
echo -en "\r$offs Head ">&2
sort -V $headers | uniq > $headers"_tmp"

egrep "^#CHROM" $headers"_tmp" > $chrompos	#Place #CHROM row at bottom, not at top
egrep "^#CHROM" -v $headers"_tmp" > $rest

cat $rest $chrompos > $headers
rm $chrompos $rest $headers"_tmp"


#Place headers + data
echo -en "\r$offs Join ">&2
cat $headers $outfile"_tmp" > $outfile
rm $headers $outfile"_tmp"
rm -rf $tmp
echo -en "\r$offs Done --> `basename $outfile` " >&2

