#!/bin/tcsh

# we pass the program arguments to the script

set seed=0
set compile=0
set gcccompile=0

while (1)
	switch ($1:q)
		case -N:
			shift
			set dim=$argv[1]
			breaksw
		case -D:
			shift
			set con=$argv[1]
			breaksw
		case -S:
			shift
			set seed=$argv[1]
			breaksw
		case -c:
			set compile=1
			breaksw
		case -C:
			set gcccompile=1
			breaksw
		case --:
			shift
			break
		default:
			echo Error!
			exit 2
		endif
	endsw
	shift
end

if ( gcccompile == 1 ) then

	set libdir=/usr/local/lib
	set incdir=/usr/local/include
	set prgnam=prg
	set cflags=(-Wall -O2 -lgsl -lgslcblas -std=c99)
	
	gcc -I $incdir -L $libdir main.c -o $prgnam $cflags

endif

# run the program ...
./prg $dim $con $seed

set mfile="matrix-$dim-$con.txt"
set ofile="connectivity-$dim-$con.txt"
set graph="graph-$dim-$con"
set mgraph="mgraph-$dim-$con"
set pfile="plot-$dim-$con.tex"

# replace the stuff in the template -- from the template we create the script
sed -e 's|CONNECTIVITY|'"$ofile"'|' \
	-e 's|NNODES|'"$dim"'|' \
	-e 's|NCONNECTIOS|'"$con"'|' \
	-e 's|MGRAPH|'"$mgraph"'.tex|' \
	-e 's|GRAPH|'"$graph"'.tex|' \
	-e 's|MATRIX|'"$mfile"'|' < plot_template.gp > $pfile

gnuplot $pfile
epstopdf $graph.eps
epstopdf $mgraph.eps

echo "$pfile created!"
echo "Plots are $graph.tex and $mgraph.tex."

# flags come into play
if ( $compile == 1 ) then
	pdflatex joze_zobec_111.tex
endif

exit 0

