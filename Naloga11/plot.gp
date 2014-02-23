set term epslatex color solid
set xlabel 'dele\v z zamenjanih povezav -- $\sigma$'

set output "plot1.tex"
set ylabel 'premer grafa -- $r$'
set title 'Premer grafa'
set log x
plot "progress-N320-D8.txt" u (2*$1/(320*8)):2 w l title '$r$'
unset output

set output "plot2.tex"
unset log
set ylabel 'gru\v cavost grafa -- $C$'
set title 'Gru\v cavost grafa'
set log y
plot "progress-N320-D8.txt" u (2*$1/(320*8)):3 w l title '$C$'
unset output

set xlabel 'indeks poteze -- $t$'
set ylabel 'dele\v z ljudi, ki so sli\v sali govorico -- $\jmath$'
set title 'Dele\v z ljudi, ki pozna govorico v odvisnosti od \v casa -- $\jmath(t)$'
unset log
set tics out
set tics nomirror
set output "joke.tex"
plot "joke-320-8-1.00.txt" u 1:2 w l title '$\sigma = 1.0$', \
"joke-320-8-0.80.txt" u 1:2 w l title '$\sigma = 0.8$', \
"joke-320-8-0.60.txt" u 1:2 w l title '$\sigma = 0.6$', \
"joke-320-8-0.40.txt" u 1:2 w l title '$\sigma = 0.4$', \
"joke-320-8-0.20.txt" u 1:2 w l lt 8 title '$\sigma = 0.2$', \
"joke-320-8-0.10.txt" u 1:2 w l lt 9 title '$\sigma = 0.1$', \
"joke-320-8-0.00.txt" u 1:2 w l lt -1 title '$\sigma = 0.0$'

unset output
