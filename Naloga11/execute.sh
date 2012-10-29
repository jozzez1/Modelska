#!/bin/tcsh

# we pass the program arguments to the script
set dim=($1)
set con=($2)
set seed=0

if ( $# == 3 ) then
	set seed=($3)
endif

# run the program ...
./prg $dim $con $seed

set mfile="matrix-$1-$2.txt"
set ofile="connectivity-$1-$2.txt"
set gfile="graph-$1-$2.tex"
set pfile="plot-$1-$2.tex"

# replace the stuff in the template
sed 's|CONNECTIVITY|'"$ofile"'|' < plot_template.gp > $pfile
sed 's|GRAPH|'"$gfile"'|' $pfile >> $pfile
sed 's|NNODES|'"$1"'|' $pfile >> $pfile
sed 's|NCONNECTIOS|'"$2"'|' $pfile >> $pfile

exit 0

