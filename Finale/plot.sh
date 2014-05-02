#!/bin/tcsh

set radii = {20,40,60,80,100}

gnuplot -p << EOF

#set term epslatex color solid size 12cm,10cm
set term pngcairo size 12cm,10cm
set out "$2.png"

set xrange [0:110]
set yrange [-20:20]
set xlabel 'oddaljenost planeta od te\\v zi\\v s\\v ca: \$\\zeta\$'
set ylabel 'radialna komponenta impulza: \$p_\\zeta\$'
set title  'Re\\v sitev na se\v cni ploskvi, \$\\p_\\psi(0) = $1\$'

plot "z${radii[1]}-L$1-N3600.dat" u 2:4 w l title '', \
	"z${radii[2]}-L$1-N3600.dat" u 2:4 w l title '', \
	"z${radii[3]}-L$1-N3600.dat" u 2:4 w l title '', \
	"z${radii[4]}-L$1-N3600.dat" u 2:4 w l title '', \
	"z${radii[5]}-L$1-N3600.dat" u 2:4 w l lt 8 title ''

unset out

EOF

exit 0
