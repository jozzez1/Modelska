#!/bin/sh

filename=`echo $1 | sed 's/-predict.dat//g'`
output=`echo $1 | sed 's/.dat//g'`
N=`wc -l $1 | awk '{print $1}'`

gnuplot << TOEND

reset

set ylabel '\$|x_t - \\hat{x}_t|\$'
set xlabel '\\v casovni indeks: \$t\$'
set xrange [0:$N]

set term epslatex color solid size 12cm,18cm
set output "$output.tex"
set key t l

set multiplot layout 2,1
set title 'Absolutna napaka linearne napovedi za \\texttt{$filename}'
plot "$1" u 1:(abs(\$2 - \$3)) w l title 'napaka pri \$p = 30\$'#, \
#	"$1" u 1:4 w l lt -1 title '\$\sigma_{\varepsilon}\$'

set title 'Primerjeva napovedi s pravimi vrednostmi'
plot "$1" u 1:2 w l title 'prave vrednosti', \
	"$1" u 1:3 w l lt 3 lw 2 title 'napoved s \$p = 30\$'

unset multiplot
unset output

TOEND

epstopdf $output.eps

exit 0

