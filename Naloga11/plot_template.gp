reset

set term epslatex color solid

set output "GRAPH"

set xlabel 'Indeks iteracije -- $i$'
set ylabel 'Premer grafa -- $d$'
set y2label 'Koeficient gru\v{c}avosti -- $c$'

plot "CONNECTIVITY" u 1:2 w l lt 1 axes x1y1, \
	"CONNECTIVITY" u 1:3 w l lt 3 axes x1y2

unset output

