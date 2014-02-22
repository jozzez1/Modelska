reset

set xlabel 'rang mesta -- $R$'
set ylabel '\v stevilo prebivalcev -- $R$'
set title 'Svetovna mesta'

f(x) = A * (log(A) / (x + log(A)));
g(x) = B * (log(B) / (x + log(B)));
h(x) = C * (log(C) / (x + log(C)));
i(x) = D * (log(D) / (x + log(D)));

fit f(x) "mestaSLO.txt" u 1:2 via A
fit g(x) "mestaZDA.txt" u 1:2 via B
fit h(x) "mestaIND.txt" u 1:2 via C
fit i(x) "mestaCHN.txt" u 1:2 via D

set log xy

#set term epslatex color solid
#set output "mest2.tex"
plot "mestaSLO.txt" lt 1 title '', "cities.txt" u 1:2 w l lt 1 title 'Slovenija', \
"mestaZDA.txt" lt 2 title '', "cities.txt" u 1:3 lt 2 w l title 'ZDA', \
"mestaIND.txt" lt 3 title '', "cities.txt" u 1:4 lt 3 w l title 'Indija', \
"mestaCHN.txt" lt 4 title '', "cities.txt" u 1:5 lt 4 w l title 'Kitajska'
#unset out

