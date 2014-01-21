# This is the configuration file for the entire pipeline.

# Touch anything apart from the triple asterisked areas ('***')
# and I will hunt you down and murder you.
#   -- Mehmet

#[Annotations IDs]
genelist_id="AL"
zygosity_id="ZYG"
type_id="TYP"
direction_id="DIL"
codon_id="COL"
protein_id="PRL"
mutation_id="MUL"
codonnom_id="CCH"
proteinnom_id="PCH"


#[Shared variables]
export delim="::::"         # Default delimiter, used numerously
export last_folder=".last"
mkdir -p $last_folder

export UNAFFECT="unaffecteds"
export AFFECT="affecteds"
export ALL="all"


#[ ***Configurable opts *** ]  <-- Go nuts.
export db="hg19"
export scenario="recessive" 			# or "dominant"
export isoforms="--isoforms"  			# (--isoforms || "" )
export vcffopts="" 						# --filter <effect> --unfilter <effect> --depth 0 -k
export bamzyopts=""  					# --qualim 0 --hetlim 60 --range 1 --extra

keepintronintergen="yes" 					# "yes"

export exgeneopts=$isoforms 			#"--exons "
export dbopts="--exons --splice 5 --scores --frames --direction" 		# --local
export commonopts=" "$isoforms  		#"--exons --isoforms"

export isointrosplops="" 				#"--no-splice"
if [ "$isoforms" = "" ];then
     export isointrosplops=$isointrosplops" --no-isoforms" 				#(too)
fi
export genetable_flags="--no-exons"

#export keywops="MISSENSE|$mutation_id|1 NONSENSE|$mutation_id|1" # &&|| SYNONYMOUS   
#                          -- deprecated by gene_tables_HTML.py, filter_funcannot.py no longer needed



#[Static data]
static_dir="/home/remote/bin/auto/z_modules/static_data"
export dnamap=$static_dir"/dna_codon.table"
export fasta_dir=$static_dir"/FASTA/$db/"
#export bam_fold=$static_dir"/BAM/"
export logfile="runlog.txt"



#[Non-configurable opts]
export dbmap="$db/$db.genemap"

export geneidopts="--geneid="$genelist_id 		#weird.. throws errors on '-', and double quoting does not work
export typeopts="--typid="$type_id				# and here. Still works though
export zygoopts="--zygid="$zygosity_id			
export zygovennopts=$zygoopts" --ignore-lowhet"

export gpopts="--multiple"
if [ "$keepintronintergen" != "" ];then
   export gpopts=$gpopts" --keepall"
fi

export matchglopts=$geneidopts
export hetchethomopts=$geneidopts" "$zygoopts" "$typeopts" "$isoforms
export funcannotopts=$dbmap" "$dnamap" "$fasta_dir" "$geneidopts" "$typeopts
export genetableopts=$geneidopts" --mutid="$mutation_id" --codid="$codonnom_id" --protid="$proteinnom_id" "$zygoopts" --genelist="$dbmap

export allopts="--ids="$genelist_id"+"\
$direction_id"+"$codon_id"+"$mutation_id"+"\
$protein_id"+"$codonnom_id"+"$proteinnom_id


#[variable denoting order of 'steps']
# The most 'efficient' pipeline is this order:
# x g f z i k n c l t 
#     Or
#
# extract.sh $vcf_fold $extr_fold
# addgenes.sh $addg_fold
# filter.sh $filt_fold $AFFECT # 3rd arg tells it what to fetch
# addzygo.sh $zygo_fold $AFFECT
# intersect.sh $inter_fold $AFFECT
# knockout.sh $knock_fold $AFFECT
# informative.sh $inf_fold
# common.sh $commfold $AFFECT
# linkage_all.sh $link_fold $AFFECT
# functional.sh $func_fold $AFFECT
# functional_genetables.sh $gtab_fold $AFFECT

