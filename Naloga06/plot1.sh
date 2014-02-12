#!/bin/sh

conf=$1

gnuplot << SCRIPT

reset

set tics nomirror
set tics out
set log y

set y2tics
set xrange [0:1024]
set key top left

set xlabel 'indeks konfiguracije -- \$i\$'
set ylabel '\$\\overline{\\chi^2}\$'
set y2label '\\v stevilo parametrov v modelski krivulji -- \$k\$'
set title 'Izra\\v cun \$\\overline{\\chi^2}\$ za vse mo\\v zne konfiguracije'

set term epslatex color solid size 15cm,9cm
set output "prelim.tex"

plot "prelim.txt" u 1:2 w l title '\$\\overline{\\chi^2}\$', "prelim.txt" u 1:3 w l lt 9 title '\$k\$' axes x1y2

unset output


SCRIPT

