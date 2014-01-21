################### File Functions #####################################

# Print number of data lines per file given as arg
fileLines(){
	files=$@
	for f in $files; do
		printCount $f
	done
#	echo "" >&2
}

# Print number of data lines in a file (not headers)
printCount(){
	if [ "$2" = "" ];then
		echo -en $1 >&2
		printright "#: $(numLines $1) "
	else
		printright "#: $(numLines $1) "
	fi
}

# Number of data lines
numLines(){
	echo $(egrep -v "^#" $1 -c)
}

# Takes a list of pedigree identifiers (seperated by '_') and a list of
# VCF files, and prints the files corresponding to the individual id
grabFiles() {
	numbers=`echo $1 | sed 's/_/ /g'`
	files=`echo -en ${@:2} | sed 's/ /\n/g'`

	for num in $numbers; do
		query="((^"$num")|(/"$num"))[_.][^/]*vcf"
		for ff in $files; do
			echo "$ff" | egrep "$query"
		done
	done
}


# /path/to/file/37_xydfas.vcf  --> 37
getNum(){
	file=$1
	echo `basename $file` | awk -F "." '{print $1}'| awk -F "_" '{print $1}'
}

# /path/to/file/37_xydfas.vcf  --> /path/to/file/37
# /path/to/file/37_xydfas.vcf  --> /path/to/37 | /path/37 | /37        with args 1, 2, 3, etc
getPrefix(){
	dname=`dirname $1`
	if [ $# = 2 ];then
		count=$2
		until [ $count = 0 ];do
			dname=`dirname $dname`
			count=$(( $count - 1 ))
		done
	fi
	echo "$dname/"$(getNum $1)
}

# Check if a folder of certain size exists , useful for determining whether
# a step has already been run or not -- resumable script
#Check if folder exists already, skip extracting -- resume
checkStepDone() {
	folder=$1
	lim=$2
	phrase=$3 #optional, helps select VCF files if multiple different types are given in the same folder
	
	if [ -d $folder ]; then
		size=`echo $(du -s $folder) | awk '{print $1}'`
		if [ $size -gt $lim ]; then
			if [ "$(ls $folder/*$phrase.vcf 2>/dev/null)" != "" ];then
				printline "Already processed $folder" $txtred '~' asd
				fileLines `ls $folder/*"$phrase".vcf`
				updateLastProcessed `ls $folder/*"$phrase".vcf`
				exit
			fi
		fi
	fi
}








##################### Last Files Processed Functions ###################

# Returns list of last processed Affected files
lastAffectedFiles() {
	ids=$(grabAllAffecteds)
#	echo $(echo $ids | sed 's/ /_/g') $(echo `ls $last_folder/*.vcf`) >&2
	grabFiles $(echo $ids | sed 's/ /_/g') $(echo `ls $last_folder/*.vcf`)
}



# Returns list of last processed Affected files
lastUnaffectedFiles() {
	ids=$(grabAllUnaffecteds)
	grabFiles $(echo $ids | sed 's/ /_/g') $(echo `ls $last_folder/*.vcf`)
}


# Returns files that match group (e.g. unaffecteds, affecteds)
# 
getGroupFiles(){
	group=$1
	annotations=""

	if   [ "$group" = "$AFFECT" ]; then annotations=$(lastAffectedFiles);
	elif [ "$group" = "$UNAFFECT" ]; then annotations=$(lastUnaffectedFiles);
	elif [ "$group" = "$ALL" ]; then annotations=$(lastAffectedFiles)" "$(lastUnaffectedFiles);
	else echo "Unknown group: $group" >&2; exit;
	fi
	
	echo $annotations
}


# Update affecteds and unaffecteds, removing only those that have
# the same prefixes as the incoming files
updateLastProcessed() {
	vcfs=$*
	current=$(ls $last_folder/*.vcf 2>/dev/null)
	
	echo -en "\r$txtpur" >&2
	
	for v in $vcfs; do
		v_pref=$(getNum $v)
	
		# If not empty, update
		if [ "$current" != "" ]; then
			for c in $current; do
				c_pref=$(getNum $c)
				
				if [ "$v_pref" = "$c_pref" ]; then
					rm $c				# remove current if incoming has same prefix
					break
				fi
			done
		fi

		cp $v $last_folder/;
		echo -en " →→→  $v_pref ">&2

	done
	
	printright "← updated"
	echo -e "$txtrst" >&2
}




##################### Pedigree Functions ###############################
# Prints out family ids in pedfile
grabFams(){
	pedfile=$1
	if [ ! -f $pedfile ];then
		echo "Cannot find pedfile" >&2
		exit
	fi

	line=`findFam.py $pedfile --affecteds --parents\
 | awk -F "[" '{print $2}'\
 | awk -F"]" '{print $1}'\
 | egrep "[0-9]" -n\
 | awk -F ":" '{print $1}'`
 
	allfams=`findFam.py $pedfile --affecteds --parents\
 | awk -F "Fam:" '{print $2}'\
 | awk '{print $1}'\
 | egrep "[0-9]"`

	fams=""
	for index in $line; do
		sh=`echo "$allfams" | head -$index | tail -1`
		fams=$fams" "$sh
	done
	echo $fams
}


grabChildren(){
	findFam.py $pedfile --affecteds --parents\
 | awk -F "[" '{print $2}'\
 | awk -F"]" '{print $1}'\
 | egrep "[0-9]"
}



#same as above until I change the dumb name 'grabChildren' to something more useful
grabAllAffecteds(){
	findFam.py $pedfile --affecteds --parents\
 | awk -F "[" '{print $2}'\
 | awk -F"]" '{print $1}'\
 | egrep "[0-9]"
}


grabAllUnaffecteds(){
	findFam.py $pedfile --unaffecteds --parents\
 | awk -F "[" '{print $2}'\
 | awk -F"]" '{print $1}'\
 | egrep "[0-9]"
}


grabUnaffSib(){
	findFam.py $pedfile --unaffecteds --parents\
 | awk -F "[" '{print $2}'\
 | awk -F"]" '{print $1}'\
 | egrep "[0-9]"
}


grabParents(){
	findFam.py $pedfile --affecteds --parents\
 | awk -F"Parents: " '{print $2}'\
 | grep -v "\-1"
}


################### Input Functions ####################################

validinput(){
	prompt=$1
	validopts=${@:2}
	
	valid="False"
	until [ "$valid" = "True" ]; do
		echo -en $prompt >&2; read var;
		for m in $validopts; do
#			echo "comparing: $var and $m" >&2
			if [ "$var" = "$m" ]; then
				valid="True"
				break;
			fi
		done
	done
	echo $var
}




################### Display Functions ##################################
export tab=8 # Default tab spacing
tabs -$tab >/dev/null 2>&1

#colours
export txtblk='\e[0;30m' # Black - Regular
export txtred='\e[0;31m' # Red
export txtgrn='\e[0;32m' # Green
export txtylw='\e[0;33m' # Yellow
export txtblu='\e[0;34m' # Blue
export txtpur='\e[0;35m' # Purple
export txtcyn='\e[0;36m' # Cyan
export txtwht='\e[0;37m' # White
export bldblk='\e[1;30m' # Black - Bold
export bldred='\e[1;31m' # Red
export bldgrn='\e[1;32m' # Green
export bldylw='\e[1;33m' # Yellow
export bldblu='\e[1;34m' # Blue
export bldpur='\e[1;35m' # Purple
export bldcyn='\e[1;36m' # Cyan
export bldwht='\e[1;37m' # White
export unkblk='\e[4;30m' # Black - Underline
export undred='\e[4;31m' # Red
export undgrn='\e[4;32m' # Green
export undylw='\e[4;33m' # Yellow
export undblu='\e[4;34m' # Blue
export undpur='\e[4;35m' # Purple
export undcyn='\e[4;36m' # Cyan
export undwht='\e[4;37m' # White
export bakblk='\e[40m'   # Black - Background
export bakred='\e[41m'   # Red
export badgrn='\e[42m'   # Green
export bakylw='\e[43m'   # Yellow
export bakblu='\e[44m'   # Blue
export bakpur='\e[45m'   # Purple
export bakcyn='\e[46m'   # Cyan
export bakwht='\e[47m'   # White
export txtrst='\e[0m'    # Text Reset

printhead(){
	cols=$2
	
	for i in `seq 1 $cols`; do echo -en $1 >&2; done
}

printright(){
	sentence=$1
	
	if [ "$width" = "" ];then
		width=$(echo `stty size` | awk '{print $2}')
	fi

	nchars=$(( `echo $sentence | wc -c` - 1))
	space=$(( $width - $nchars ))
	numtabs=$(( $space / $tab ))
	extraspace=$(( $space - $numtabs ))

	build="\r"
	for i in `seq 2 $numtabs`; do build=$build"\t"; done
	for j in `seq 1 $extraspace`; do build=$build" "; done

	echo -e $build$sentence >&2
}

printline(){
	color=$2
	dim=$3
	width=$COLUMNS


	if [ "$width" = "" ];then
		width=$(echo `stty size` | awk '{print $2}')
	fi
	

	if [ $# = 1 ];then
		color=$undgrn
		dim=" "
	elif [ $# = 2 ]; then
		dim=" "
	fi
	

	echo -en $color >&2
	printhead " " $width
	

	nchars=$(( `echo $1 | wc -m` + 5 ))
	lendim=$(( `echo $dim | wc -m` -1 ))
	
	if [ $nchars = 0 ]; then nchars=1;fi
	if [ $lendim = 0 ]; then lendim=1;fi
	
	offset=$((  (( $width - nchars ) / 2 ) / $lendim ))
	
	SPACE=""
	for i in `seq 1 $offset`; do SPACE=$SPACE$dim; done
	
	echo -e " $SPACE$1$SPACE " >&2
	echo $1 >> $logfile
	
	#reset
	if [ $# = 4 ];then
		echo -en $txtrst >&2
	else
		echo -e $txtrst >&2
	fi
}

#################################################################
# Linkage function
# Search for linkage (messner) file (if any)
findLinkageFile() {
	files=$(ls ./*.bed 2>/dev/null)
	linkage=""
		
	if [ "$files" = "" ]; then
		printright "No linkage bed file found.. continuing."
		linkage="no"
	else
		linkage=""
		for file in $files; do
			# multiple check
			if [ "$linkage" != "" ];then
				echo -e ":multiple BED files detected! Only one can be used for linkage:\n $files" >&2
				echo 255
				exit -1
			fi
			
			# Check if file is in <whatever>_hg<ver>.bed format --> version part is important
			version=`echo $file | egrep "hg[12][01289]" | awk -F "hg" '{print $NF}' | awk -F ".bed" '{print $1}'`
			if [ "$version" = "" ];then
				echo -e "Linkage: $file -- if this is to be used to filter for linkage, please rename it to \"<whatever>_hg<version>.bed\" \n  e.g.  messner_hg19.bed" >&2
				echo 255
				exit -1

			elif [ "$version" != "19" ]; then
				echo -en "Linkage: $file -- Linkage data is hg$version and sequencing data is most likely hg19.\nWould you like to transpose linkage positions from $version --> 19?" >&2
				ans=$(validinput "(y/n)" "y n yes no")
				if [[ "$ans" =~ "y" ]]; then
					echo -en "Converting..." >&2
					tmp_link="linkage_hg19.bed"
					messner_transform $file hg19 > $tmp_link
					linkage=$tmp_link
					echo "..made $tmp_link" >&2
				else
					linkage=$file
				fi
			else linkage=$file;
			fi
		done
	fi
	echo $linkage
}


### MISC ####
folder_count=/var/tmp/folder_count.txt				# Stores what order of processing we are on, independent of folder names
rm $folder_count;									# Reset every run.
echo "0" > $folder_count

# takes $name $number, prints 01_folder, 02_etc, etc
foldername(){
	name=$1
	count=$(cat $folder_count)
	newcount=$(( $count + 1 ))
	echo $newcount > $folder_count
	echo $(printf "%02d" $newcount)"_$name"
}

# Exports are only visible to SUB-processes, not parents. Calling source needs no exports
export -f checkStepDone
export -f grabFiles
export -f numLines
export -f fileLines
export -f printCount
export -f getNum
export -f getPrefix
export -f lastAffectedFiles
export -f lastUnaffectedFiles
export -f getGroupFiles
export -f updateLastProcessed
export -f grabFams
export -f grabChildren
export -f grabAllAffecteds
export -f grabAllUnaffecteds
export -f grabUnaffSib
export -f grabParents
export -f validinput
export -f printline
export -f printhead
export -f printright
export -f findLinkageFile
export -f foldername
