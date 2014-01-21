#!/bin/bash

usage(){
	echo "`basename $0` <VCF fold> <BAM fold>" >&2
	echo ""
	echo "where the VCF files are in gz archives and contain both snvs and indels (RAW)"
	exit
}

if [ $# != 2 ]; then
   usage
fi

vcf_fold=$1
export bam_fold=$2


################ Load Ext Modules and Make GeneMap #####################
# Find pedfile and linkage file
. initialise_config.sh

# Search for linkage file
export linkage=$(findLinkageFile)
[ "$linkage" = "255" ] && exit


######################## G E N T L E M E N ! ############################
################ S T A R T ## Y O U R ## E N G I N E S #################
clear
printline "Next-Generation Sequencing Pipeline" $txtpur " " none
echo
printline "Royal Free Hospital, Nephrology Dept.  © 2013" $txtpur " " none
printline "Bionformatics: M Tekman, M Mozere, H Stanescu" $txtpur " " none

echo ""

########################################################################
########### B E G I N  ### T H E ### P R O C E S S I N G ###############
. load_genemap.sh 


# Extract
#    Combine snvs and indels into one file for each file in the VCF folder
printline " Extract: Extracting and combining SNVs and INDELs "
$(extract.sh $vcf_fold `foldername $extr_fold`)

# Add Genes
#    Annotate the files with genes
printline " Add Genes: Annotating the VCF files"
$(addgenes.sh `foldername $addg_fold`)
#
#


# Add Zygo
#    Annotate the files with zygosity
printline " Add Zygo: Appending Zygosity to filtered VCF files "
$(addzygo.sh `foldername $zygo_fold` $AFFECT)
#
#

# Functional
#   Add functional annotation
#   Filter by these functional annotations  (if filter settings set)
#   Print Gene tables (if table settings set)
printline " Functional: Add functional annotations "
$(functional.sh `foldername $func_fold` $AFFECT)
#
#


# Functional Filter
# Remove isoforms or introns
#printline " Functional Filter: Remove isoforms or introns "
#$(functional_filter.sh `foldername $func_fold"_filter"` $AFFECT)


# Common
#   Finally produce a common genelist, and match all affecteds for the remaining variants within those genes
#printline " Common: Common genes "
#$(common.sh $commfold $AFFECT)
#
#

# Filter
#    Filter the annotated files
#printline " Filter: Filtering the Affecteds "
#$(filter.sh `foldername $filt_fold` $AFFECT) # 3rd arg tells it what to fetch
#
#


# Linkage
#   Filter by Linkage
printline " Linkage: Applying Linkage "
$(linkage_all.sh `foldername $link_fold` $linkage $AFFECT)
#
#


# Intersect
#    Matches lines for all individuals sharing common genes generated from
#    the variant overlap between siblings and the gene overlap for nonsibs
printline " Intersect: Intersecting Sib overlaps, with NonSib overlaps "
$(intersect.sh `foldername $inter_fold` $AFFECT)
#
#


# Knockout
#    From the remaining genes/lines, determing what is knockout (KO) and
#    what is non-knockout (NKO)
printline " Knockout: Determining Knockout (KO) and Nonknockout (NKO) genes "
$(knockout.sh `foldername $knock_fold` $AFFECT)
#
#


# Informative
#    Grab common variants from parents for affected individuals
#       Kick out HOM variants that are HOM for unaffected siblings too (Recessive)
#    OR Kick out any variants that are HET for unaffected siblings (Dominant)
printline " Informative: Remove uninformative variants ($scenario) "
$(informative.sh `foldername $inf_fold` )
#
#


# Common
printline " Common "
$(common.sh `foldername $comm_fold` $AFFECT)
#
#

# Gene Tables
printline " GeneTables "
$(genetables.sh `foldername $gtab_fold` $AFFECT)
#
#
