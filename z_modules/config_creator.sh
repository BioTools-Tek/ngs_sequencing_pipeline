#!/bin/bash

#Defaults
dbopts="--local --exons --splice 5"
gpopts="--multiple"
vcffopts=""
linkage="no"
bamzyopts=""

usage(){
	echo "Interactive configuration file creator for the ngs_pipeline.sh script" >&2
	echo
	echo "usage: `basename $0` [--default]"
	echo
	echo "--default    gives default config"
	exit
}

askyesno(){
	read -p "$1 [y/n]?" -n 1 ans
	if [ "$ans" = "y" ]; then
		echo $2	# true
	else
		echo $3 # false
	fi	
}

interact(){
	#Reset all
	dbopts=""; gpopts=""; vcffopts=""; linkage="no"; bamzyopts="";
	
	echo -e "Bedtarget options:\n"
	#Q1
	read -p "Use local database [y/n]?" -n 1 ans
	if [ "$ans" = "y" ]; then
		dbopts=$dbopts" --local"
	fi
	
	#Q2
	read -p "Use genes instead of exons?[y/n]?" -n 1 ans
	if [ "$ans" = "y" ]; then
		dbopts=$dbopts" --genes"
	else
		dbopts=$dbopts" --exons"
	fi
	
	#Q3
		
	
}



if [ $# -gt 0 ]; then
	if [ "$1" == "--default" ];then
		exit 0
	else
		usage
	fi

else
	interact
fi

