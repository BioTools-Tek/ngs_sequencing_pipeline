#!/bin/bash

usage(){
	echo "
	where files are seperated by $delim
	
	`basename $0` <[VCF_FILE_LIST]> <venn_output> <output_fold>
OR
	`basename $0` <[VCF_FILE_LIST]> <venn_output> <output_fold> <venn_accepts> <venn_rejects>
	    if you want a zygosity check too
	    (make sure to load configuration script first)
	" >&2
	exit
}

if [ $# -lt 3 ];then
	usage
fi


files=`echo $1 | sed "s/$delim/ /g"`
venn_output=$2
numfold=$3
venn_accepts=$4
venn_rejects=$5

mkdir -p $numfold

venn.py $files --columns --overlap --sep "$delim" > $venn_output

if [ $# = 5 ];then
	echo "Filtering out questionable zygosities" >&2
	zygosity_venn_check.py $venn_output $venn_rejects --delim=$delim $zygovennopts > $venn_accepts
else
	venn_accepts=$venn_output
fi

venn2vcf.sh $venn_accepts --delim="$delim$delim" --folder=$numfold		#twice delim
echo `ls $numfold/*venned*`
