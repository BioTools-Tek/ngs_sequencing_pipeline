#!/bin/bash

# 2013_07_30 -- changed venn combine file, to just combine combined file

if [ "$delim" = "" ];then		# . load_config.sh
	delim="::::"
fi

p_arg="--parents="

usage(){
	echo "`basename $0` <output fold> [$p_arg<file1.vcf>[+<file2.vcf>]] [Affected Child VCFs ONLY]>
		
		at least one parent must be given. Two parents must be seperated by a \"$delim\",
			e.g. --parents=mother.vcf$delimfather.vcf
			
		Children must be listed (with spaces) at the end
	
" >&2
	exit
}

if [ $# -lt 3 ];then
	usage
fi

fold=$1
vcfs=${@:2}
parents=""

if [[ "$2" =~ "$p_arg" ]]; then
	parents=`echo $2 | sed "s/$p_arg//" | sed "s/$delim/ /"`
	vcfs=${@:3}
fi

num_chil=`echo $vcfs | wc -w`

if [ "$num_chil" = "0" ]; then
	echo "No Affected children present" >&2
	exit
fi

mkdir -p $fold
if [ "$parents" = "" ]; then
	if [ "$num_chil" = "1" ]; then
		name=$(echo `basename $vcfs` | sed 's/.vcf/_venned.vcf/')
		echo "$vcfs -- only 1 file given, nothing to venn, copying over to $fold/$name" >&2
		cp $vcfs $fold/$name
		exit
	fi
fi

########################################################################

# Venn the affecteds against each parent
venned_list=""   # affecteds may appear more than once depending on num parents
for parent in $parents; do
	numpar=$(getNum $parent)
	parfold=$fold/"parents"/$numpar"_and_affchildren"/
	overlap="$parfold/parent_affchild.overlap"

	venned=$(venn_and_split.sh `echo "$vcfs $parent" | sed "s/ /$delim/g"` $overlap $parfold)
	rm $overlap
	
	#Grab affecteds only
	parentfile=$(grabFiles $numpar $venned)
	affectfile=`echo -e $venned | sed "s/ /\n/g" | grep -v $parentfile`
	
	#copy parents as dummies
	newname=$(getPrefix $parentfile 3)"_parentdummy.vcf"
	cp $parent $newname
	
	venned_list=$venned_list" "$affectfile
done


# From child-mother and child-father, produce unique variant list for child
# (i.e. do not venn, sort variants)
for vcf in $venned_list; do
	num=$(getNum $vcf)
	files=$(grabFiles $num $venned_list)
	
	ifold=$fold/$num"_combined"/

	#If we have 2 parents, then 2 files to merge per affected
	if [ "`echo $files | wc -w`" = "2" ];then
		mkdir -p $ifold
		tmpname=$ifold/$(basename $vcf)
		newname=$(getPrefix $vcf 3)"_inherited.vcf"
				
		combinefiles.sh $tmpname $files;
		echo "" >&2
		cp $tmpname $newname
	else
		#For single parent affected childen, they've already been combined
		newname=$(getPrefix $vcf 3)"_inherited.vcf"
		cp $vcf $newname
	fi
	#remove the VCF file from the venned list
	echo $newname
	venned_list=`echo -e $venned_list | sed "s/ /\n/g" | grep -v $vcf`
done
