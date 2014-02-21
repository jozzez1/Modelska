reset

set xlabel 'indeks to\v cke -- $N$'
set ylabel 'amplituda signala'

set term epslatex color solid

set output "res-signal0.tex"
set title '{\tt signal0.dat}'
plot "signal0.dat" w l title '$c(t)$', \
"results0.dat" u 1:2 w l lt -1 title '$u(t)$'
unset output

set output "res-signal1.tex"
set title '{\tt signal1.dat}'
plot "signal1.dat" w l title '$c(t)$', \
"results1.dat" u 1:2 w l lt -1 title '$u(t)$'
unset output

set output "res-signal2.tex"
set title '{\tt signal2.dat}'
plot "signal2.dat" w l title '$c(t)$', \
"results2.dat" u 1:2 w l lt -1 title '$u(t)$'
unset output

set output "res-signal3.tex"
set title '{\tt signal3.dat}'
plot "signal3.dat" w l title '$c(t)$', \
"results3.dat" u 1:2 w l lt -1 title '$u(t)$'
unset output

