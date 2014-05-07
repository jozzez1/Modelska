#!/bin/tcsh

gnuplot -p << EOF

set term epslatex color solid size 12cm,10cm
set out "$1.tex"

set xlabel 'polarni kot planeta: \$\\psi\$'
set ylabel '\$\\zeta\$ pri \$p_\\psi = 160\$'
set y2label '\$\\zeta\$ pri \$p_\\psi = 220\$'
set tics nomirror
set tics out
set y2tics

set key t l

set xrange [0:2*pi]
set xtics ("0" 0, '\$\\pi/2\$' 0.5*pi, '\$\\pi\$' pi, '\$3\\pi/2\$' 1.5*pi, '\$2\\pi\$' 2*pi)

plot "stable_orbit.txt" u 5:4 w l title '\$p_\\psi = 160\$', \
	"stable_orbit2.txt" u 5:4 w l axes x1y2 title '\$p_\\psi = 220\$' lt 3

unset out

EOF

exit 0

