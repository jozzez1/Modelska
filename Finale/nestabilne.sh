#!/bin/tcsh

gnuplot -p << EOF

set term epslatex color solid size 9cm,9cm
set out "$2.tex"

set polar
set grid polar

plot "$1" u 5:4 w l title 'planet', "$1" u 3:(\$2/4000) w l title '\$M_1\$', "$1" u 3:((-1)*\$2/2000) w l title '\$M_2\$'

EOF

exit 0
