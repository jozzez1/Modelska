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
set efile="graph-$1-$2.eps"
set mgraph="mgraph-$1-$2.tex"
set mefile="mgraph-$1-$2.eps"

# replace the stuff in the template -- from the template we create the script
sed 's|CONNECTIVITY|'"$ofile"'|' < plot_template.gp > tmpfile
sed 's|GRAPH|'"$gfile"'|' < tmpfile > $pfile
sed 's|NNODES|'"$1"'|' < $pfile > tmpfile
sed 's|NCONNECTIOS|'"$2"'|' < tmpfile > $pfil
sed 's|MGRAPH|'"$mgraph"'|' < $pfil > tmpfile
sed 's|MATRIX|'"$mfile"'|' < tmpfile > $pfil

rm tmpfile

gnuplot $pfile

epstopdf $efile
epstopdf $mefile

exit 0

