#!/bin/tcsh

gnuplot -p << EOF

set term epslatex color solid size 12cm,10cm
set out "$1.tex"

set key b r

f(x) = k*x
g(x) = K*x

set xlabel 'kvadrirana vrtilna koli\\v cina planeta: \$p_\\psi^2\$'
set ylabel 'zastojni radij: \$\\zeta_0\$'
set title  'Obna\\v sanje zastojnih to\\v ck, \$\\varepsilon = 0\$'

fit f(x) 'equilibrium_e0.txt' u (\$1**2):2 via k
fit g(x) 'equilibrium2-M4000_e0.txt' u (\$1**2):2 via K

plot 'equilibrium_e0.txt' u (\$1**2):2 w p title '\$M_1 = 2000\$', \
	f(x) lt -1, \
	'equilibrium2-M4000_e0.txt' u (\$1**2):2 w p title '\$M_1 = 4000\$', \
	g(x) lt 3

unset out

EOF

exit 0

