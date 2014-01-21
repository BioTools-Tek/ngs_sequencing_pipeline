#!/bin/bash

# Valid command keys:
declare -A commands=( 
	#actual names					#shortcuts
	["extract"]="extract"			["e"]="extract"
	["filter"]="filter"				["f"]="filter"
	["addgenes"]="addgenes"			["g"]="addgenes"
	["addzygo"]="addzygo"			["z"]="addzygo"
	["intersect"]="intersect"		["i"]="intersect"
	["knockout"]="knockout"			["k"]="knockout"
	["informative"]="informative"	["n"]="informative"
	["common"]="common"				["c"]="common"
	["linkage"]="linkage"			["l"]="linkage"
)
# Valid argument keys:
declare -A arguments=( 
	#actual names					#shortcuts
	["affected"]="affected"			["a"]="affected"
	["unaffected"]="unaffected"		["u"]="unaffected"
)


usage(){
	echo "Takes in a list of sequencing commands (and their suffix arguments) and executes them sequentially in the order given, where the output from each process is piped into the next
	
	usage: `basename $0` <[command strings]>
	
	type `basename $0` --help     for a list of valid keywords and arguments
	"; exit
}

help(){
	build=""
	for keys in "${!commands[@]} ${!arguments[@]}"; do
		build=$build" "$keys
	done
	echo -e $build
#	echo $build | sort | uniq
}

help
