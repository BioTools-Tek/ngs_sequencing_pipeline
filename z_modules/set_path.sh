#!/bin/bash

usage(){
	echo "Exports scripts as system wide without duplicating PATH variables each time
	
	usage: `basename $0` <location of auto folder>
	"
	return 0
}

if [ $# != 1 ]; then
	usage
fi

autoscripts=$1
savepath="$HOME/.savedpath";

if [ -f $savepath ];then
	#Load old path
	PATH=`head -1 $savepath`
else
	echo $PATH > $savepath
fi

# Find autoscripts, make executable

scripts=`find -L $autoscripts -type f | grep -v "~" | sort`
for script in $scripts; do chmod a+x $script; done

# Find autoscript directories, add to path
autonames=$(find -L $autoscripts -type d | grep -v ".git")
export PATH=$PATH`echo $autonames | sed 's/ /:/g'`
