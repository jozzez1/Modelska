#!/bin/tcsh

gnuplot -p << EOF

set term epslatex color solid
set output "pentlja.tex"

set xrange [-10:210]
set yrange [-40:40]

set xlabel 'kartezi\\v cna koordinata \$x\$'
set ylabel 'kartezi\\v cna koordinata \$y\$'

set polar
set grid polar

plot "stable_orbit3.txt" u 5:4 w l title ''

unset out

EOF

exit 0
