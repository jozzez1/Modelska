reset
set term epslatex color solid
set out 'clovestvo.tex'
set xlabel 'leto'
set ylabel 'svetovna populacija [$10^6$]'
set key bot
set tics out
set tics nomirror
set title '\v Clove\v ska populacija po Kapitsi'
plot "zgodovina1.dat" u 1:2 w p title '\texttt{zgodovina.dat}', \
"zgod1.fit" u 1:2 w l lt -1 title 'fit od 1800', \
"zgod0.fit" u 1:2 w l lt 4 title 'celoten fit'

unset out
