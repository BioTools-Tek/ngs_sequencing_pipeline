#!/bin/bash


usage(){
	echo "`basename $0` <output folder>"
	
	exit
}

if [ $# != 1 ];then
	usage
fi

folder=$1
aff_files=$(getGroupFiles $AFFECT)

#This WONT work for dominant scenario (parents will be affected too)
unaff_children=$(grabUnaffSib $pedfile)
unaff_parents=$(echo `grabParents $pedfile` | sed 's/[+\n]/_/g')

unaffecteds=`echo $unaff_children" "$unaff_parents | sed s'/\s/_/g'`
unaffected_vcfs=$(lastUnaffectedFiles)



# Processing each entire unaffecteds file takes a looong time.
# Instead extract a genelist (be as general as possible) from the affected
# and match gene lines with the unaffecteds, and then combine based upon that subsselection of useful variants.
matchfold=$folder/"a_narrow_unaffecteds"
un_arg=$(echo $unaffected_vcfs | sed "s/ /$delim/g")
af_arg=$(echo $aff_files | sed "s/ /$delim/g")
unaffected_vcfs=$(narrow_unaffecteds.sh $matchfold $af_arg $un_arg)

allfiles=$aff_files" "$unaffected_vcfs

# Combine affecteds with their parents to get inherited variants.
printline " Finding inherited variants " $txtylw "~"
filtered1_all=$(handle_parents_all.sh $pedfile $folder/"b_family_variants" $allfiles)

# Combine affecteds with their unaffected siblings to kick out useless variants
printline " Combining affecteds and unaffected siblings " $txtylw "~"
filtered2=$(handle_sibs_all.sh $pedfile $bam_fold $folder/"c_"$scenario"_sibling_variant_filter" $filtered1_all)

#echo $filtered2
updateLastProcessed $filtered2



