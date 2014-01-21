#!/bin/bash

usage(){
	echo -e"
`basename $0` <output folder>
	
 Performs genepender annotations on all files in .last/ folder
 Batch caller for genepender.
 
 Please initialise config before running."
	exit -1
}

if [ $# -lt 1 ]; then
	usage
fi

folder=$1
vcfs=$(ls $last_folder/*.vcf)

rejfold="rejects"

#checkStepDone $folder 200   # Individual file check now in progress

processThese(){
	annotations="";
	fold=$1
	rej=$fold/$2
	files=${@:3}
	
	if [ "$files" = "" ];then
		echo "Nothing to process" >&2
	else
		mkdir -p $rej
		
		for vcf in $files; do
			rejects=$rej/$(getNum $vcf)"_unwanted.vcf"
			appender=$fold/$(getNum $vcf)"_annotation.vcf"
			if [ -f $appender ];then
				echo -en "$appender " >&2
			else
				echo $appender >&2
				genepender $dbmap $vcf $rejects $gpopts > $appender
			fi
			printCount $appender "no"
			annotations=$annotations" "$appender
		done
		echo $annotations
	fi
}

affected_ids=`echo $(grabChildren) | sed 's/[ \n]/_/g'`
affecteds=$(grabFiles $affected_ids $vcfs)

# Only works for RECESSIVE scenario where parents are unaffected
unaffected_ids=`echo $(grabUnaffSib) | sed 's/[ \n]/_/g'`
unaffected_ids=$unaffected_ids"_"`echo $(grabParents) | sed 's/[ \n+]/_/g'`
unaffecteds=$(grabFiles $unaffected_ids $vcfs)

printline " Processing Affecteds " $txtblu "~" asd
affs=$(processThese $folder/$AFFECT $rejfold $affecteds)

printline " Processing Unaffecteds " $txtblu "~" asd
unaffs=$(processThese $folder/$UNAFFECT $rejfold $unaffecteds)

updateLastProcessed $affs $unaffs >&2
