#!/bin/tcsh

# we pass the program arguments to the script

set seed=0
set compile=0
set gcccompile=0
set dim=0
set con=0
set help=0

while ($#argv > 0)
	switch ($argv[1])
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
		case -latex:
			set compile=1
			echo Setting up the latex compile.
			breaksw
		case -gcc:
			set gcccompile=1
			echo Setting up the compilation of prg.
			breaksw
		case -h:
			set help=1
			breaksw
		case --:
			shift
			breaksw
		default:
			shift
			break
		endif
	endsw
	shift
end

if ( $dim == 0 || $con == 0 ) then
	if ( $help == 0 ) then
		echo Arguments N and D needed!
		echo Correct usage:
		echo "$0 -N <number of nodes> -D <number of connections>"

		exit 1
	else
		echo
		echo "Obligatory flags for $0"
		echo "-------------------------------"
		echo "-N <n>"
		echo "	n is number of nodes you want for your graph"
		echo
		echo "-D <d>"
		echo "	d is number of connections that each node has"
		echo "	in the beginning -- must be an even number"
		echo
		echo "Optional flags for $0"
		echo "-------------------------------"
		echo "-S <s>"
		echo "	s is the seed value of the Marsenne Twister"
		echo "	pseudo-random generator."
		echo 
		echo "-gcc"
		echo "	compile the source main.c and output it in the"
		echo "	program prg"
		echo
		echo "-latex"
		echo "	compile the report joze_zobec_111.tex with"
		echo "	pdflatex"
		echo
		echo "-h"
		echo "	print this message"

		exit 0
	endif
endif

if ( $gcccompile == 1 ) then

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

