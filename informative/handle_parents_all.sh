#!/bin/bash

usage(){
	
	echo "`basename $0` <pedfile.pro> <output folder> <[ALL VCFs]>"
	echo ""
	exit
}

if [ $# -lt 3 ];then
	usage
fi

pedfile=$1
folder=$2
allfiles=${@:3}

mkdir -p $folder

#. load_functions.sh
fams=$(grabFams $pedfile)
affchildren=$(grabChildren)
unaffchildren=$(grabUnaffSib)
parents=$(grabParents)

#num_fams=`echo $fams | wc -w`

fam_filtered=""

fam_num=0
for fam_id in $fams; do
	fam_num=$(( $fam_num + 1 ))
	logfam=$folder/$fam_id"_log"
	#
	affs=`echo $affchildren | awk -v ind=$fam_num {'print $ind'}`
	unaffs=`echo $unaffchildren | awk -v ind=$fam_num {'print $ind'} | sed s'/_/ /g'`
	pars=`echo $parents | awk -v ind=$fam_num '{print $ind}' | sed s'/+/_/g'`
	#
	par_vcfs=`echo $(grabFiles $pars $allfiles) | sed "s/ /$delim/g"`								#<parent1>::::<parent2>
	aff_vcfs=$(grabFiles $affs $allfiles)								#<Aff1>  <Aff2>
	#
	unaff_ids=`echo $unaffs | sed s'/ /_/g'`
	unaff_vcfs=$(grabFiles $unaff_ids $allfiles)

	echo -e "Family: $fam_id,  Affs[$affs]   \tParents[$pars]" >&2
	
	for u in $unaff_vcfs; do
		cp $u $folder/$(getNum $unaff_vcfs)"_sibdummy.vcf";
	done

	mendelian_suspects=""
	if [ "$par_vcfs" = "" ];then # if no parents, then just copy the affecteds
		echo "       No parent VCFs found" >&2
		for affchild in $aff_vcfs; do
			newname=$folder/$(getNum $affchild)"_affected.vcf"
			cp $affchild $newname
			mendelian_suspects=$mendelian_suspects" "$newname
		done
	else
		fold=$folder/"Fam_"`echo $fams | awk -v ind=$fam_num '{print $ind}'`

		# Copy over unaffecteds into output folder (will be used in handle_sibs_all.sh step
		mkdir -p $fold
		
	#	mendelian_suspects=$(overlap_parents_with_affecteds.sh "family_variants_$fold" $family_vcfs 2>/dev/null)
		mendelian_suspects=$(overlap_parents_with_affecteds.sh $fold --parents=$par_vcfs $aff_vcfs)
	fi
#	echo $mendelian_suspects >&2
	fam_filtered=$fam_filtered" "$mendelian_suspects
done

echo $fam_filtered
