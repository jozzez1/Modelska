reset

set term epslatex color solid

set output "GRAPH"

set title '$N = NNODES$, $D = NCONNECTIOS$'
set tics nomirror

set yrange [0:5.2]
set y2range [0:0.5]

set ytics 1
set y2tics 0.1

set xlabel 'Indeks iteracije -- $i$'
set ylabel 'Premer grafa -- $d$'
set y2label 'Koeficient gru\v{c}avosti -- $c$'

plot "CONNECTIVITY" u 1:2 w l lt 1 axes x1y1 title '$d$', \
	"CONNECTIVITY" u 1:3 w l lt 3 axes x1y2 title '$c$'

unset output

