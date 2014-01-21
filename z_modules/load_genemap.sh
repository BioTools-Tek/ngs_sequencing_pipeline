#Supported flags
maps="hg18 hg19"

if [ "$dbmap" = "" ] || [ ! -f "$dbmap" ];then
	printline " GeneMap not found, creating new one " $txtred "~" asd
	db=$(validinput "Human Genome version [$maps]: " $maps)

	printline "Creating GeneMap $db"

	if [ "$dbopts" = "" ];then
		echo "initialise config please " >&2
		exit
	fi
	makegenemaps.sh $db $dbopts

else
	printline " GeneMap found $dbmap " $txtblu "~" asd
fi
