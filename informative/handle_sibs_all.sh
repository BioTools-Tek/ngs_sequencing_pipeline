#!/bin/bash

usage(){
	echo "
	`basename $0` <pedfile.pro> <BAM folder> <output fold> <[All VCF files in pedigree]>
	
	" >&2
	exit -1
}

if [ $# -lt 4 ];then
	usage
fi

pedfile=$1
bam_fold=$2
folder=$3
fams=$(grabFams $pedfile)
files=${@:4}

#	For each affected output VCF from previous step, kick out any
#	homozygous variants that the affected VCF file shares with its
#	unaffected siblings

grabFromFam(){
	fam_id=$1
	affect=$2
	query="Fam:$fam_id[[:space:]]$affect"
	res=$(findFam.py $pedfile | egrep "$query")
	echo `echo $res | sed "s/$query//" | awk -F "Sibs:" '{print $2}'`
}

filtered=""


for fam_id in $fams; do
	logfam=$folder/$fam_id"_log"

	affsib=$(grabFromFam $fam_id "Aff")
	unaffsib=$(grabFromFam $fam_id "Unaff")
	
	affnums=`echo $affsib | sed 's/_/ /g'`
	unaffnums=`echo $unaffsib | sed 's/_/ /g'`
	
	for aff in $affnums; do
		aff_file=$(grabFiles $aff $files)
		if [ "$aff_file" = "" ];then
			echo "No VCF file for $aff" >&2
		else
			echo -en "Family: $fam_id,\tComparing $aff  ">&2

			outfold=$folder/$aff"_filter"
			mkdir -p $outfold
			
			if [ "$unaffnums" = "" ];then
				echo -en "--> No sibs!" >&2
				cp $aff_file $outfold/
				aff_file=$outfold/`basename $aff_file`
			else
				for unf in $unaffnums; do
					unf_file=$(grabFiles $unf $files)
					if [ "$unf_file" = "" ];then
						echo -en "--> ($unf missing) " >&2
						cp $aff_file $outfold/
						aff_file=$outfold/`basename $aff_file`
					else
						# Each filter through an unaffected sib should work on the output from the previous filtering of the previous sibling.
						# so the affected file should be updated from each sib
						echo -en "--> $unf " >&2
#						aff_file=$(knockout_sib_variants.sh $aff_file $unf_file $bam_fold $outfold 2>>$logfam)
						aff_file=$(knockout_sib_variants.sh $aff_file $unf_file $bam_fold $outfold)
					fi
				done
			fi
			echo "--> `basename $aff_file`" >&2
		fi
		filtered=$filtered" "$aff_file
	done
done

echo $filtered
