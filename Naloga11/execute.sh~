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
sed -e 's|CONNECTIVITY|'"$ofile"'|' \
	-e 's|NNODES|'"$1"'|' \
	-e 's|NCONNECTIOS|'"$2"'|' \
	-e 's|MGRAPH|'"$mgraph"'|' \
	-e 's|GRAPH|'"$gfile"'|' \
	-e 's|MATRIX|'"$mfile"'|' < plot_template.gp > $pfile

gnuplot $pfile
epstopdf $efile
epstopdf $mefile

echo "$pfile created!"
echo "Plots are $efile and $mefile."

# flags come into play

exit 0

