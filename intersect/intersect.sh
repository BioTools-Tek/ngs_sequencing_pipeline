#!/bin/bash

NONE="NONE"

usage(){
	echo "`basename $0` <output_folder> <group>
	
 Matches lines for all individuals sharing common genes generated from
 the variant overlap between siblings and the gene overlap for nonsibs
"
	exit -1
}

if [ $# != 2 ]; then
	usage
fi

folder=$1
group=$2

zygosity_files=$(getGroupFiles $group)

checkStepDone $folder 200 matchedlines

common="overlapping_genes"

s4sib=$folder/$common/"sibling"
s4nonsib=$folder/$common/"nonsibling"
common_genes=$folder/$common/"common_genes.txt"

# Step 4 
# a. Perform variant comparisons between affected siblings and grep combined genes
common_sib_genes=$(sibling_genelists.sh $s4sib $pedfile $zygosity_files)
#
# b. Grep combined genes of affected non-sibs
# Provide a red ferrari.
####### T H I S IS A CLONE OF SIBLING_GENELISTS AND DOES NOT MAKE SENSE ########
common_nonsib_genes=$(nonsibling_genelists.sh $s4nonsib $pedfile $zygosity_files)

#
# c. Overlap all genes from above
printline " Overlapping.." $txtpur "_" asd

# Check if any data from previous step
if [ "$common_sib_genes" = "$NONE" ] && [ "$common_nonsib_genes" = "$NONE" ]; then
	echo "No data -- Problem" >&2
	exit -1
fi

if [ "$common_sib_genes" = "$NONE" ] || [ "$common_nonsib_genes" = "$NONE" ]; then
	if [ "$common_sib_genes" = "$NONE" ]; then
		cat $common_nonsib_genes > $common_genes
	fi
	
	if [ "$common_nonsib_genes" = "$NONE" ]; then
		cat $common_sib_genes > $common_genes
	fi
else
	overlap_genes.py $common_sib_genes $common_nonsib_genes > $common_genes
fi

#
# d. Match lines for those genes for all files
echo -e "\nExtracting genelines.." >&2
matched=""
for vcf in $zygosity_files; do
	newname=$folder/$(echo `basename $vcf` | sed 's/.vcf/_matchedlines.vcf/' )
	match_gene_lines.py $vcf $common_genes $matchglopts > $newname
#	echo $newname >&2
	matched=$matched" "$newname
done

#echo $matched
updateLastProcessed $matched
