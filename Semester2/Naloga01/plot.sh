#!/bin/sh

filename=`echo $1 | sed 's/-spectra.dat//g'`
output=`echo $1 | sed 's/.dat//g'`

gnuplot << TOEND

reset
#set log y

set ylabel '\$\log_{10}|P(\w)|^2\$'
set xlabel 'frekvenca: \$\w = 2\pi\nu/N\$'

set xtics ('0' 0, '\$\pi/2\$' pi/2, '\$\pi\$' pi, '\$3\pi/2\$' 3*pi/2, '\$2\pi\$' 2*pi);
set xrange [0:2*pi];
set key b r

set term epslatex color solid size 12cm,16cm
set output "$output.tex"

set multiplot layout 2,1 title '\\texttt{$filename}'
plot "$1" u 1:(log10(\$2)) w l title '\$p = 10\$', \
	"$1" u 1:(log10(\$3)) w l title '\$p = 20\$', \
	"$1" u 1:(log10(\$4)) w l title '\$p = 30\$', \
	"$1" u 1:(log10(\$5)) w l title '\$p = 40\$'

plot "$1" u 1:(log10(\$6)) w l title '\\texttt{$filename}', \
	"$1" u 1:(log10(\$4)) w l lt 3 lw 2 title '\$p = 30\$'

unset multiplot
unset output

TOEND

epstopdf $output.eps

exit 0

