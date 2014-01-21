#!/bin/bash

# Caller for kickuninformative.py

usage(){
	echo "Removes variants from affected VCF that also exist in unaffected VCF with unhelpful zygosities
	
	`basename $0` <affected.vcf> <unaffected.vcf> <BAM fold> <output fold> (recessive|dominant)
	
	" >&2
	exit -1
}

if [ $# -lt 4 ]; then
	usage
fi

aff=$1
unaff=$2
bam_fold=$3
folder=$4

if [ "$5" != "" ];then
	if [ "$5" = "recessive" ]; then
		scenario="recessive"
	elif [ "$5" = "dominant" ];then
		scenario="dominant"
	else
		echo "Please specify either 'dominant' or 'recessive'" >&2
		exit
	fi
fi

mkdir -p $folder

uid=$( echo `basename $unaff` | awk -F"_" '{print $1}')
aid=$( echo `basename $aff` | awk -F"_" '{print $1}')

# Venn siblings together, and extract only the unaffected
# Apply bamzygo to subselection of unaffected data
venn_folder=$folder/$aid"_venn"/
venn_overlap=$venn_folder/"venn.overlap"
bamzy_fold=$folder/$aid"_bamzygo"/

overlapped_files=$(venn_and_split.sh $aff"$delim"$unaff $venn_overlap $venn_folder)
unaffected_subsample=$(grabFiles $uid $overlapped_files)

echo "BAMMING UP" $unaffected_subsample >&2
unaffected=$(addzygo.sh $bamzy_fold --bamfold=$bam_fold $unaffected_subsample)


#zygovennopts="--zygid=ZYG"   # . load_config.sh
#scenario="recessive"        # . load_config.sh

newname=$folder/$(echo `basename $aff` | sed s"/.vcf/_sibout$uid.vcf/")

kicksibvariants.py $aff $unaffected $zygovennopts $scenario > $newname
mv `ls $aid*rejects\.vcf` $folder/$rejname

echo $newname
