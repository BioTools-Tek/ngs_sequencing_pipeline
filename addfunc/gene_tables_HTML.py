#!/usr/bin/env python

import sys
from VCF_Operations import VCFOps
from time import ctime, sleep

version="2013_09_19"

#args
g_arg="--geneid="
m_arg="--mutid="
c_arg="--codid="
p_arg="--protid="
z_arg="--zygid="
ne_arg="--no-exons"
gl_arg="--genelist="

def usage():
	print >> sys.stderr, '''
v%s
Takes in multiple functionally annotated VCF files and prints an HTML table showing the mutations that occured for each gene for each patient

	usage: %s %sAL %sMUL %sCCH %sPCH %sZYG %s<genelist> [OPTIONS] [<VCF Files>]

OPTIONS:
%s	groups genes by gene name and not exons
''' % (version,\
sys.argv[0].split('/')[-1], g_arg, m_arg, c_arg, p_arg, z_arg,\
gl_arg,\
ne_arg)
	exit(-1)



## Parse Args ##
if len(sys.argv)<8:
	usage()

g_id=sys.argv[1].split(g_arg)[1]
m_id=sys.argv[2].split(m_arg)[1]
c_id=sys.argv[3].split(c_arg)[1]
p_id=sys.argv[4].split(p_arg)[1]
z_id=sys.argv[5].split(z_arg)[1]
genelist=sys.argv[6].split(gl_arg)[1]

if ( "" in (g_id, m_id, c_id, p_id, z_id ) ):
	print >> sys.stderr, "Require Ids in the order given:\n"
	usage();

noexons=False
files=[]

for arg in sys.argv[7:]:
	if arg.startswith(ne_arg):
		noexons=True
	else:
		files.append(arg)


class GeneContainer:
	def __init__(self, genename, chrom, pos1, pos2):
		self.name = genename;
		self.chrom = chrom
		self.pos1 = pos1
		self.pos2 = pos2



########################################################################
print >> sys.stderr, "Making Gene Map:",
gene_positions={}

numlines = VCFOps.countFileLines(genelist)
linecount=0

f = open(genelist,'r')
#chr1	69090	69090	OR4F5|Exon1_SpliceD	cC	+	0
for line in f:
#	print >> sys.stderr, line
	linecount += 1
	
	tokens = line.split('\t')
	
	chrom = tokens[0].strip()
#	print >> sys.stderr, tokens
	pos1 = int(tokens[1])
	pos2= int(tokens[2])
	gene=tokens[3].strip()
	
	if (gene.find("Splice")!=-1) or (gene.find("Intron")!=-1):
		continue
	
	if noexons:
		gene=gene.split('|')[0].strip();
	
	if gene in gene_positions:
		gc = gene_positions[gene]
		
		for gg in gc:
			if (gg.chrom == chrom):
				gg.pos1 = pos1 if (pos1 < gg.pos1 ) else gg.pos1
				gg.pos2 = pos2 if (pos2 > gg.pos2 ) else gg.pos2
			else:
				gene_positions[gene].append(GeneContainer(gene,chrom,pos1,pos2))
#		else:
#			print >> sys.stderr, "Disagreement in gene:", gene
#			exit(-1)
	else:
		gene_positions[gene] = [GeneContainer(gene, chrom, pos1, pos2)]
#		print >> sys.stderr, "Added:", gene
	
	if linecount%90==0:
		print >> sys.stderr, "\r\t\t", (100*linecount)/numlines,'%',

print >> sys.stderr, "\r\t\t100 %     "
#print >> sys.stderr, [(x.chrom, x.pos1, x.pos2,) for x in gene_positions["AARS2"]]
#exit(-1)


########################################################################
print >> sys.stderr, "Making Mutation Map:",
########################################################################
def getType(codon_change):
	if (codon_change.find("ins")!=-1):
		return "ins"
	if (codon_change.find("del")!=-1):
		return "del"
	return "snv"

numtabs=2;
gene_map={}      # Populated across all files
                 # map -> Gene -> File -> List of codon and protein changes
chrom_map={}     # Gene-->Chrom

arg_ids = [g_id, m_id, c_id, p_id]   #z_id omitted because it doesn't have the same number of elements

for file in files:
	
	file=file.strip()
	prefix=file.split('/')[-1].split('_')[0]

	vcf = VCFOps(file)
	vcf.getFieldIndexes(arg_ids+[z_id])
	
	f = open(file,'r')
	
	countlines=0

	for line in f:
		countlines += 1
		
		if line[0]=='#':
			continue
		
		func_data = line.split('\t')[vcf.IFORMAT_INDEX].split(':')
		data_map = vcf.getElementsFromData(func_data, check=arg_ids);
		
		chrom = line[0:5]
		zyg_data = func_data[vcf.field_map[z_id]].split(',');
		zygosity = zyg_data[-1][0:3];
		read_depth = int(zyg_data[-2][1:])
		
		for index in xrange(len(data_map[g_id])):
			
			gene = data_map[g_id][index]
			mutation = data_map[m_id][index]
			codon_change = data_map[c_id][index]
			protein_change = data_map[p_id][index]
			
			#uniq = (   (c. 14A>G, p. 3C>R) , ( snv, MIS, HET, 40 ) )
			uniq = ( (codon_change, protein_change), ( getType(codon_change), mutation[0:3], zygosity, read_depth ) )

			if noexons:
				gene = gene.split('|')[0]

			if gene in gene_map:
				if file in gene_map[gene]:
					gene_map[gene][file].append( uniq )
				else:
					gene_map[gene][file] = [ uniq ]
			else:
				gene_map[gene] = {}
				gene_map[gene][file] = [ uniq ]
				
			
			if gene not in chrom_map:
				chrom_map[gene]=chrom;

		if (countlines%21==0):
			print >> sys.stderr, '\r', '\t'*numtabs, prefix,':', int(100*countlines/vcf.numlines), '%',
	
	print >> sys.stderr, '\r', '\t'*numtabs, "%s:100%%   " % prefix,
	f.close()
	numtabs += 1





########################################################################
print >> sys.stderr, "\nPrinting HTML:",
########################################################################

genecount=0;                        # progress
mapkeys=sorted(gene_map.keys());
numgenes = len(mapkeys)     # num genes

tablename="genetable"

#classes
mutie_class="mutie"
type_class="types"
data_class="data"
syn_class="syn"
patient_class="pat"
common_class="common"
depth_class="depth"
chrom_class="chrom"
genename_class="genename"

#ids -- this may perhaps be pointless...
toggmute_but = "toggmute_but"
toggtype_but = "toggtype_but"
toggsyns_but = "toggsyns_but"
toggblanks_but = "toggblanks_but"
toggsort_but = "toggsort_but"
depthabove_but = "depthabove_but"

#styles
disc_style="font-size:0.6em"

#HTML HEADERs
print '''
<meta charset="UTF-8">
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">'''

print '''<!-- \nThe following is autogenerated by gene_tables_HTML.py\nc. Mehmet Tekman @ Royal Free Hospital, Nephrology\n --!>'''
print "<head>"


## CSS
print '''
<style>
#buttons { font-size:0.7em; }
#buttons input {
	border: 2px solid;
	border-radius: 5px;
	width:100px;
	height:25%;
}

#'''+tablename+''' {
	font-family:courier, "courier new", monospace;
	border-collapse:collapse;
}

#'''+tablename+''' td, #'''+tablename+''' th {
	word-wrap: break-word;
	font-size:0.8em;
	border:1px solid #98bf21;
	padding:3px 7px 2px 7px;
}

#'''+tablename+''' th {
	font-family:"Trebuchet MS", Arial, Helvetica, sans-serif;
	font-size:90%;
	text-align:center;
	padding-top:5px;
	padding-bottom:4px;
	background-color:#A7C942;
	color:#fff;
}

#genecounter {
	position:fixed;
	font-size:1em;
	bottom:0px;
	left:0px;
	background-color:#000000;
	color:white;
}

.gene {
	color:#111111; 
	font-weight: bold;
	text-align:center;
	width:100px;
}

.'''+data_class+''' .'''+type_class+''' {color:blue;}
.'''+data_class+''' .'''+mutie_class+''' {display:inline-block;}
.'''+depth_class+''' {color:purple; width:30px; border:1px solid #ffccff; border-radius:5px; float:right; display:inline-block;}
.'''+chrom_class+''' {font-size:0.8em; color:#FF6666; font-family:courier, "courier new", monospace;}

.nasil { 
	display:inline-block;
	border: 1px solid #dddddd;
	border-radius:3px 3px 3px 3px;
	width:150px;
}

.common {background-color:#ccddcc; max-width:200px;}
.common_genes {
	border-bottom:1px solid #ddeedd;
	border-right:1px solid #ddeedd;
//background-color:#aaffaa;
	border-radius:10px;
}



</style>
'''

#JAVASCRIPT
print '''
<script>
//Global Pointers
var tabla="", gcounter="", num_affs=-1, depth_box="";
var default_butt_color='#fff0ff', inactive_butt_color='#999999';
var button_dict = {}, class_dict = {};


function initialise(){
	tabla = document.getElementById("'''+tablename+'''");
	gcounter = document.getElementById("genecounter");
	depth_box = document.getElementById("'''+depthabove_but+'''");

	depth_box.style = "width:30px;"
	depth_box.value = "0";
	
	class_dict = {"'''+toggmute_but+'''":"'''+mutie_class+'''","'''+toggtype_but+'''":"'''+type_class+'''","'''+toggsyns_but+'''":"'''+syn_class+'''"}
	button_dict = {"'''+toggmute_but+'''":true, "'''+toggtype_but+'''":true, "'''+toggsyns_but+'''":true, "'''+toggblanks_but+'''":false, "'''+toggsort_but+'''":false}
	
	//Set default colours
	for (k in button_dict) setButtonColor(k)
	
	num_affs = tabla.rows[0].cells.length-2;
	checkBlanks_UpdateCommon();
}

// ------- Hiding attributes for the nasil container object --------- //

//The following flags define which classes independently affect visibility of nasil
//If >1 are active in the object, then nasil is hidden. If ALL inactive, nasil is shown
var nasil_flags= { "'''+syn_class+'''":'synflag', "'''+depth_class+'''":'depthflag' } // extensible

function initFlag(object_ref){
	for (k in nasil_flags)
		if (typeof object_ref[nasil_flags[k]] === "undefined") 
			object_ref[nasil_flags[k]]=false;
}

function flag_toggle(object_ref,booler, flag_class){
	initFlag(object_ref)
	object_ref[nasil_flags[flag_class]] = booler
}
function valueFlag(object_ref,flag_class){
	initFlag(object_ref)
	return object_ref[nasil_flags[flag_class]];
}
//untested
function allflagsDown(object_ref){
	var all_false = true;
	for (k in nasil_flags){
		if(object_ref[nasil_flags[k]]){
			all_false = false;
			break;
		}
	}
	return all_false;
}

function flagDEPTHOn (object_ref){ flag_toggle(object_ref, true , "'''+depth_class+'''");}
function flagDEPTHOff(object_ref){ flag_toggle(object_ref, false, "'''+depth_class+'''");}

//Unused
//function flagSYNOn (object_ref){ flag_toggle(object_ref, true , "'''+syn_class+'''");}
//function flagSYNOff(object_ref){ flag_toggle(object_ref, false, "'''+syn_class+'''");}


// ---- Button Functions ---------------------------------------------//
/** Logs button presses and handles colors **/
function buttonPress(button_id){
	var booler = (button_dict[button_id] = !button_dict[button_id]);
	setButtonColor(button_id);
	return booler
}

/** Switches display for specific class objects, or changes the 
 *  hide value of the nasil container object **/
function toggleClass(button_id, incrementparent=false, checkblanks=false, checkpatients=true)
{
	var booler = buttonPress(button_id),  class_name = class_dict[button_id];

	var cells=tabla.getElementsByClassName(class_name)
	for (var s=0; s < cells.length; s++){
		var parent = cells[s];

		if (incrementparent){
			while( parent.className != "nasil" ) parent = parent.parentNode;
			flag_toggle(parent, !booler, class_name);
		} 
		else parent.style.display = booler?"":"none";
	}
	if (checkpatients) togglePatients();
	if (checkblanks) checkBlanks_UpdateCommon()
}


function togglePatients(){
	var hide_bool = ( (!button_dict["'''+toggmute_but+'''"]) && (!button_dict["'''+toggtype_but+'''"]) )
	cells=tabla.getElementsByClassName("'''+data_class+'''")
	for (var c=0; c< cells.length; c++) cells[c].style.display=hide_bool?"none":"";
	return hide_bool;
}

function setButtonColor(button_id){
	document.getElementById(button_id).style.background = button_dict[button_id]? default_butt_color : inactive_butt_color;
}

function toggleSyn(){ toggleClass("'''+toggsyns_but+'''", incrementparent=true, checkblanks=true, checkpatients=false) }
function toggleMut(){ toggleClass("'''+toggmute_but+'''") };
function toggleTyp(){ toggleClass("'''+toggtype_but+'''") };
function toggleBlank(){
	buttonPress("'''+toggblanks_but+'''")
	checkBlanks_UpdateCommon()
};


// Filtering function upon read depth

function setDepth(){
	curr_depth = +depth_box.value;  //unary plus treats as numeric
	cells=tabla.getElementsByClassName("'''+depth_class+'''")

	for (var s=0; s < cells.length; s++){
		var parent = cells[s], innervalue = parent.innerHTML.trim();
		while( parent.className != "nasil" ) parent = parent.parentNode;
		
		if ( innervalue <= curr_depth ) flagDEPTHOn(parent) //incrementHide(parent)
		else flagDEPTHOff(parent)
	}
	checkBlanks_UpdateCommon()
}


/**  1. Updates row (tr element) visibility if all 'nasil' mutation containers are hidden
 *   2. Find common mutations across rows
 *   3. returns number of visible gene rows 
 **/
function checkBlanks_UpdateCommon(){
	var num_genes=0;
	var blank_value = button_dict["'''+toggblanks_but+'''"]?num_affs:1;
	//ShowBlanks ? true = { show row iff #empties >= 1 } : false = { show row iff #empties >= #num_affs }

	for (var y=1, row; row = tabla.rows[y]; y++)
	{
		var cells = row.cells, last = cells.length-1;
		var idname=cells[0].innerHTML.trim();
		
		var repeated_vals = {};
		var empty_cells_in_row=0

		for (var x=1, cell; cell = cells[x]; x++)
		{
			if (x==last) {
				var outkeys=[]
				for (k in repeated_vals) if (repeated_vals[k] >= num_affs) outkeys.push("<span class='common_genes' >"+k+"</span>")
				cell.innerHTML = (outkeys.length==0)?"-":outkeys.join("<br/>");
				break
			}
			else {
				var cell_empty=true;
				var spanners = cell.getElementsByClassName("nasil"); //each mutation
				
				for (var s=0; s < spanners.length; s++)
				{
					var nasil = spanners[s];
				
					mute_obj = nasil.getElementsByClassName("'''+mutie_class+'''")[0];
					type_obj = nasil.getElementsByClassName("'''+type_class+'''")[0];
					
					// If both of the main containers are hiding, increment hide nasil
					if ( (mute_obj.style.display=="none") && (type_obj.style.display=="none") ) incrementHide(nasil);

					//If nothing is set, then show nasil
					if ((!valueFlag(nasil,"'''+syn_class+'''")) && (!valueFlag(nasil,"'''+depth_class+'''"))){
						var mutie_name = mute_obj.innerHTML.trim();
					
						//Add displayed to common pool
						if (!(mutie_name in repeated_vals)) repeated_vals[mutie_name] = 1;
						else repeated_vals[mutie_name] +=1;
						
						cell_empty=false;
						nasil.style.display="";
					}
					else nasil.style.display="none";
				}
				if (cell_empty) empty_cells_in_row ++;
			}
		}
		if (empty_cells_in_row >= blank_value ) row.style.display="none";
		else {
			row.style.display=""
			num_genes ++;
		}
	}
	gcounter.innerHTML = "&nbsp"+num_genes+" genes&nbsp";
}


// --- Sorting Table by Chromosome, GeneName -------------------------//
function sortTable(){
	var sort_chrom = buttonPress("'''+toggsort_but+'''")
	
    var store = [];
    var rows = tabla.rows
    for(var i=1, len=rows.length; i<len; i++){
        var row = rows[i];
        var sortnr = null;
        
        if (sort_chrom) sortnr = parseInt(row.cells[0].getElementsByClassName("'''+chrom_class+'''")[0].innerHTML.split("chr")[1]);
		else { // sort gene_name
			sortnr = row.cells[0].getElementsByClassName("'''+genename_class+'''")[0].innerHTML;
		}
		store.push([sortnr, row]);
    }
    if (sort_chrom){
		store.sort(function(x,y){
			return x[0] - y[0];
		});
	} else store.sort()
	
    for(var i=1, len=store.length; i<len; i++){
        tabla.appendChild(store[i][1]);
    }
    store = null;
}


</script>'''

print "</head>"


#instructions
print "<body onload=\"initialise()\" >"
print "<span style=\""+disc_style+"\">Toggle data using buttons below</span>"


#Genecounter (fixed)
print "<div id=\"genecounter\"></div>"

### Begin Table ###
print "<table id=\""+tablename+"\">"

#First row
print "<tr><th>"


## Buttons
print '''
<div id="buttons">
<input type="button" id="'''+toggmute_but+'''" value="Nomencl." onclick="toggleMut()" /><br />
<input type="button" id="'''+toggtype_but+'''" value="Detail" onclick="toggleTyp()" /><br />
<input type="button" id="'''+toggsyns_but+'''" value="Synonymous" onclick="toggleSyn()" /><br/>
<input type="button" id="'''+toggblanks_but+'''" value="Show Blanks" onclick="toggleBlank()" /><br/>
<input type="button" id="'''+toggsort_but+'''" value="Sort Chrom" onclick="sortTable()" /><br/>
Depth Above: <input type="number" id="'''+depthabove_but+'''" value="0" onchange="setDepth()" />
</div>'''
print "</th>"

# patient ID
for file in files:
	print "<th class=\""+data_class+"\" align=center>",file.split('/')[-1],"</th>"
print "<th class=\""+data_class+" "+common_class+"\" style=\"font-size:0.8em;text-align:center;\" >Common<br/>Mutations</th></tr>"


#Subsequent rows
elm_sep=" "  #tab;
mut_sep="<br/>"


for gene in mapkeys:
	gene=gene.strip()
	
	print "<tr id=\"",gene,"\" >"
	print "<td class='gene' ><span class=\""+genename_class+"\" >",gene,"</span><br />"
	
	#Find the gene that is in the same chromosome
	chr_gene = chrom_map[gene].strip()
	gene_name_dupes = gene_positions[gene]
	
	gene_pos=""
	for g in gene_name_dupes:
		if g.chrom == chr_gene:
			gene_pos = g.chrom+":"+str(g.pos1)+"<br />-"+str(g.pos2)
			break;  #found one in the same chrom, good enough or sorting will be a mess..
	
	print "<span class=\""+chrom_class+"\" >"+gene_pos+"</span></td>"
	
	for file in files:
		print "<td class=\"",data_class,"\" >"
		try:
			for key in gene_map[gene][file]:
				cp, detail = key
				
				c,p = cp
				snv,typ,zyg,depth = detail
				
				nomen=c+elm_sep+p
				next_space_index = nomen.index(' ');

				# if nomenclature is long, split at every 20 characters
				maxlen=18
				if next_space_index>maxlen:
					tmp_nom=""
					while len(nomen)>maxlen:
						tmp_nom += nomen[0:maxlen]+" "
						nomen = nomen[maxlen:]
					nomen = tmp_nom;
				
				misc_classes = " "
				if typ=="SYN":
					misc_classes += syn_class+" "
				
				print "<div class=\"nasil\" >"
				print "<span class=\"",mutie_class,"\" >",nomen,"</span>"
				print "<br />"
				print "<span class=\"",type_class+"\" >",\
							"<span ",(("class=\""+misc_classes+"\" >") if (len(misc_classes)>3) else "\" >"),\
								snv,elm_sep,typ,elm_sep,zyg,elm_sep,\
							"</span>",\
							"<span class=\"",depth_class,"\">",\
								str(depth),\
							"</span>",\
						"</span>"
				print "</div>"
				print mut_sep

		except KeyError:
			print " " # file doesnt have that gene
		
		print "</td>"
	print "<td class=\"",data_class," ",common_class,"\" align=center ></td>"

	genecount += 1
	print "</tr>"
	
	print >> sys.stderr, '\r\t\t', int(100*genecount/numgenes), '%',
print >> sys.stderr, "\r\t\t100%  "
print "</table>"

#Disclaimer
print "<span style=\""+disc_style+"\">"
print "GeneTables_v"+version+" - NGS Pipeline (<i>Nephrology, Royal Free Hospital</i>)<br/>"
print "<b>  Authors: </b>Horia Stanescu, Mehmet Tekman, Monika Mozere @ <i>University College London</i><br/>"
print "<b>Generated: </b> "+ctime()+"<br/><br/>"
print "note: Gene positions are derived directly from the gene map, which may be configured to print coding regions only."
print "</span>"

print "<br />"*3
print "</body>"
print "</html>"
