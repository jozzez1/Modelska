#!/bin/sh

gnuplot << TOEND

reset

set term epslatex color solid size 12cm,10cm
set xlabel 'kiralnost -- \$x\$'
set ylabel 'dol\v zina disklinacij -- \$\\ell\$'

set output "solve-eight.tex"
set title '{\\tt wrboost\\_theta.dat}'
plot "wrboost_eight.dat" u 1:2:(2) w yerror title "eksperiment" , "609-1.txt" w l lt 3 title 'regresijska krivulja'
unset output

set output "solve-omega.tex"
set title '{\\tt wrboost\\_omega.dat}'
plot "wrboost_omega.dat" u 1:2:(2) w yerror title "eksperiment", "522-2.txt" w l lt 3 title 'regresijska krivulja'
unset output

set output "solve-theta.tex"
set title '{\\tt wrboost\\_theta.dat}'
plot "wrboost_theta.dat" u 1:2:(2) w yerror title "eksperiment", "536-3.txt" w l lt 3 title 'regresijska krivulja'
unset output

TOEND

exit 0

