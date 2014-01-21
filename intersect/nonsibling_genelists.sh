#!/bin/bash

# Produces a combined genelist of unrelated indivuals
#

usage(){
	echo "`basename $0` <output folder> <pedfile> <[VCF files]>"
	exit -1
}

if [ $# -lt 3 ]; then
	usage
fi

folder=$1
pedfile=$2
files="${@:3}"

nonsib=$folder/"affected_nonsibs.txt"
comgen=$folder/"common_genes_affected_nonsiblings.txt"

mkdir -p $folder
printline "Step 4b: NonSibling specific...." $txtpur "_" asd

findFam.py $pedfile --affecteds --nonsiblings > $nonsib
nonsibs=`cat $nonsib | awk -F "Sibs:" '{print $2}'`

if [ "`echo $nonsibs | wc | awk '{print $2}'`" = "0" ]; then
	echo "No 'non siblings'" >&2
	exit 0
fi


gene_list=""

for nums in $nonsibs; do
	nums=`echo $nums | sed 's/[\n+ ]/_/'`

	numfold=$folder/$nums/
	mkdir -p $numfold

	#Find vcf files corresponding to each affected individual
	grep_files=$(grabFiles $nums $files)
	grep_files=`echo $grep_files | sed "s/ /$delim/g"`
	
	venn_output=$numfold/"venn_"$nums".output"
	venn_overlap=$numfold/"venn_"$nums".overlap"
	venn_rejects=$venn_overlap".badzygosities"
	venn_accepts=$venn_overlap".acceptable"
	combined_genes=$numfold/"overlapping_genes.txt"

	#Only 1 file, nothing to venn
	if [ "`echo $nums | sed s'/_/ /g' | wc -w`" = "1" ];then
		echo "$nums has no siblings" >&2
		genes=`echo $file | sed 's/.vcf/.genes/'`
		extract_genes_vcf.py $geneidopts $grep_files $exgeneopts > $combined_genes        #combined because there's only file to combine

	else
		venned_files=$(venn_and_split.sh $grep_files $venn_output $numfold $venn_accepts $venn_rejects)

		# Extract genes for each overlapped file (should be the same, but just in case)
		genelist=""
		for file in $venned_files; do
			genes=`echo $file | sed 's/.vcf/.genes/'`
			extract_genes_vcf.py $geneidopts $file $exgeneopts > $genes
			genelist=$genelist" "$genes
		done
		overlap_genes.py $genelist > $combined_genes

	fi
	
	gene_list=$gene_list" "$combined_genes
	
done

#Take the list of sibling genelists and overlap them
overlap_genes.py $gene_list > $comgen

echo $comgen

